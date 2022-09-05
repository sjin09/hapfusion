import os
import json
import time
import math
import pysam
import bisect
import natsort
import numpy as np
import hapsmash.util
import hapsmash.cslib
import hapsmash.bamlib
import hapsmash.haplib
import hapsmash.vcflib
import multiprocessing as mp
from collections import defaultdict
from dataclasses import dataclass, field
from collections import defaultdict, Counter
from typing import Set, Dict, List, Tuple


@dataclass
class METRICS:
    ccs: int = 0
    lq_ccs: int = 0 
    hq_ccs: int = 0
    phaseable_hq_ccs: int = 0
    unphaseable_hq_ccs: int = 0
    phased_hq_ccs: int = 0 
    unphased_hq_ccs: int = 0 
    hap_consistent_hq_ccs: int = 0 
    hap_inconsistent_hq_ccs: int = 0 
    recombinant_ccs: int = 0 


def get_recombination(
    chrom: str,
    chrom_len: int,
    bam_file: str,
    vcf_file: str,
    loci_lst: List[Tuple[str, int, int]],
    min_mapq: int,
    min_alignment_proportion: float,
    min_sequence_identity: float,
    qlen_mean: float,
    qlen_lower_limit: float,
    qlen_upper_limit: float,
    min_bq: int,
    min_trim: float,
    md_threshold: int,
    mismatch_window: int,
    max_mismatch_count: int,
    chrom2recombination_lst: Dict[str, List[List[Tuple[str, int, str, str]]]],
    chrom2recombination_statistics: Dict[str, List[int]]
) -> List[Tuple[str, int, str, str, int, int, int, float, float]]:

    # counter = 0
    m = METRICS()
    alignments = pysam.AlignmentFile(bam_file, "rb")
    if os.path.exists(vcf_file) and vcf_file.endswith(".vcf"):
        hetset, homset, hetpos_set, hompos_set, pos2allele = hapsmash.vcflib.load_snp_indel(chrom, vcf_file)

    (
        phased_hpos_lst,
        phased_hetsnp_lst,
        phased_hetsnp_set,
        hblock_lst,
        hidx2bidx,
        hidx2hstate,
        hidx2hetsnp,
        hetsnp2bidx,
        hetsnp2hidx,
        hetsnp2hstate,
    ) = hapsmash.vcflib.get_phased_hetsnps(vcf_file, chrom, chrom_len)
    del hidx2bidx
    del hblock_lst
    del hidx2hstate
    del hetsnp2hidx
    del hidx2hetsnp

    hapsmash_lst = []
    # hapsmash_ccs_lst = []
    # hapsmash2coord = defaultdict(list)
    # hetsnp2cnt = defaultdict(lambda: 0)
    seen = set()
    for loci in loci_lst:
        chunkloci_lst = hapsmash.util.chunkloci(loci)
        for chunkloci in chunkloci_lst:
            chunk_start, chunk_end = chunkloci[1:]
            allelecounts = defaultdict(lambda: np.zeros(6))
            if vcf_file.endswith(".bgz"):
                hetset, homset, hetpos_set, hompos_set, pos2allele = hapsmash.vcflib.load_bgz_snp_indel((chrom, chunk_start - qlen_mean, chunk_end + qlen_mean), vcf_file)
           
            for i in alignments.fetch(*chunkloci):
            # for line in alignments.fetch(*chunkloci):
                m.ccs += 1
                ccs = hapsmash.bamlib.BAM(i)
                hapsmash.util.update_allelecounts(ccs, allelecounts)
                if not ccs.is_primary: 
                    m.lq_ccs += 1
                    continue

                if ccs.mapq < min_mapq:
                    m.lq_ccs += 1
                    continue

                if ccs.qlen < qlen_lower_limit or ccs.qlen > qlen_upper_limit:
                    m.lq_ccs += 1
                    continue
               
                if ccs.query_alignment_proportion < min_alignment_proportion:
                    m.lq_ccs += 1
                    continue
                
                if hapsmash.util.get_blast_sequence_identity(ccs) < min_sequence_identity:
                    m.lq_ccs += 1
                    continue
                
                m.hq_ccs += 1
                ccs.load_mutations(hetset, homset, hetpos_set, hompos_set, phased_hetsnp_set, pos2allele) 
                if len(ccs.hetsnp_lst) == 0: # homozygous region 
                    m.unphaseable_hq_ccs += 1 
                    continue 
                elif len(ccs.hetsnp_lst) == 1: 
                    m.phaseable_hq_ccs += 1
                    continue
                else: 
                    m.phaseable_hq_ccs += 1 
                    
                ccs_hbit_lst, h0_hbit_lst, h1_hbit_lst, phased_hetsnp_subset_lst = hapsmash.haplib.get_ccs_hbit_lst(
                    ccs,
                    hetsnp2bidx,
                    hetsnp2hstate,
                    phased_hpos_lst,
                    phased_hetsnp_lst,
                )
                if len(ccs_hbit_lst) == 0: # ccs belongs to more than one haplotype block
                    m.unphased_hq_ccs += 1 
                    continue 
                m.phased_hq_ccs += 1 
                if ccs_hbit_lst == h0_hbit_lst: 
                    m.hap_consistent_hq_ccs += 1    
                    continue
                elif ccs_hbit_lst == h1_hbit_lst: 
                    m.hap_consistent_hq_ccs += 1    
                    continue 

                m.hap_inconsistent_hq_ccs += 1
                ccs.hap = hapsmash.haplib.get_ccs_haplotype(h0_hbit_lst, ccs_hbit_lst) 
                ccs_hapsmash_hetsnp_candidate_lst = hapsmash.haplib.get_hapsmash_hetsnps(ccs.hap, h0_hbit_lst, ccs_hbit_lst, phased_hetsnp_subset_lst) # search candidate 
                if len(ccs_hapsmash_hetsnp_candidate_lst) == 0: ## hetsnp is deleted ## is this the result of MMR or sequencing error?
                    continue
                   
                ccs_hapsmash_hetsnp_lst = []
                trimmed_qstart = math.floor(min_trim * ccs.qlen)
                trimmed_qend = math.ceil((1 - min_trim) * ccs.qlen)
                for (tpos, ref, alt) in ccs_hapsmash_hetsnp_candidate_lst:
                    
                    qpos = ccs.tpos2qpos[tpos]
                    _, qbq = ccs.tpos2qbase[tpos]
                    if qbq < min_bq:
                        continue

                    if qpos < trimmed_qstart:
                        continue
                    if qpos > trimmed_qend:
                        continue

                    mismatch_start, mismatch_end = hapsmash.util.get_mismatch_range(tpos, qpos, ccs.qlen, mismatch_window)
                    jdx = bisect.bisect_left(ccs.mismatch_tpos_lst, mismatch_start)
                    kdx = bisect.bisect_right(ccs.mismatch_tpos_lst, mismatch_end)
                    mismatch_count = kdx - jdx 
                    if mismatch_count > max_mismatch_count:
                        continue

                    ins_count = allelecounts[tpos][4]
                    del_count = allelecounts[tpos][5]
                    total_count = sum(allelecounts[tpos])
                    if del_count != 0 or ins_count != 0:
                        continue
                    
                    if total_count > md_threshold:
                        continue
                    
                    ccs_hapsmash_hetsnp_lst.append((tpos, ref, alt))

                if ccs_hapsmash_hetsnp_lst != 0:
                    if ccs.qname in seen:
                        continue
                    
                    seen.add(ccs.qname)
                    tend_lst = []
                    tstart_lst = []
                    hapsmash_ccs_lst = [] 
                    hapsmash_region_lst = []
                    print("{}:{}-{}".format(ccs.tname, ccs.tstart, ccs.tend)) 
                    for j in alignments.fetch(ccs.tname, ccs.tstart, ccs.tend):
                        read = hapsmash.bamlib.BAM(j)
                        # print(read.qname)
                        if read.is_primary and read.mapq >= min_mapq: 
                            hapsmash.cslib.cs2tpos2qbase(read)
                            read.load_mutations(hetset, homset, hetpos_set, hompos_set, phased_hetsnp_set, pos2allele) 
                            read_hbit_lst, h0_hbit_lst, h1_hbit_lst, phased_hetsnp_subset_lst = hapsmash.haplib.get_ccs_hbit_lst(
                                read,
                                hetsnp2bidx,
                                hetsnp2hstate,
                                phased_hpos_lst,
                                phased_hetsnp_lst,
                            )
                            # if len(phased_hetsnp_lst)
                            read.hap = hapsmash.haplib.get_ccs_haplotype(h0_hbit_lst, read_hbit_lst) 
                            if len(read_hbit_lst) == 0: # read belongs to more than one haplotype block
                                continue 
                            if read_hbit_lst == h0_hbit_lst: 
                                continue
                            elif read_hbit_lst == h1_hbit_lst: 
                                continue 
                            read_hapsmash_hetsnp_candidate_lst = hapsmash.haplib.get_hapswitch_hetsnps(read.hap, h0_hbit_lst, read_hbit_lst, phased_hetsnp_subset_lst) 
                            

                            
                            print("{}:{}\t{}:{}-{}".format(read.qname, read.hap, read.tname, read.tstart, read.tend))
                            print(h0_hbit_lst)
                            print(h1_hbit_lst)
                            print(read_hbit_lst)
                            print(phased_hetsnp_subset_lst)
                            print() 
                        # print(read.qname, read.is_primary, read.mapq, read.qlen, qlen_lower_limit, qlen_upper_limit)
                        # tend_lst.append(ccs.tend)
                        # tstart_lst.append(ccs.tstart)
                    # hapsmash_ccs_lst.append((ccs.qname, ccs.tstart, ccs.tend))
                        # tend_lst.append(ccs.tend)

                    # hetsnp2cnt[tpos, ref, alt] += 1
                    # hapsmash_ccs[tpos, ref, alt].append([ccs.tstart, ccs.tend])
                    # ccs_hapswitch_hetsnp_lst.append((chrom, tpos, ref, alt))
   
    # for (ccs.qname, ccs.tstart, ccs.tend) in hapsmash_ccs_lst:
    #     if ccs.qname in seen:
    #         continue
        ## select all CCS reads 
        ## determine start and end of the diagram
        ## determine the mutations 
        ## determine CCS reads with and without gene conversion or crossover
        ## determine ccs_hbit for each ccs read 
        ## determine consensus h0_hbit
        ## determine cosnensus h1_hbit
         
        
        # coord_lst = hapsmash2coord[(tpos, ref, alt)]
        # if len(coord_lst) == 1:
        #     start = coord_lst[0][0]
        #     end = coord_lst[0][1]
        # else:
        #     start_lst = []
        #     end_lst = []
        #     for coord in coord_lst:
        #         start_lst.append(coord[0])
        #         end_lst.append(coord[1])
        #     start = min(start_lst)
        #     end = max(end_lst)

            
        # for line in alignments.fetch(*(chrom, start, end)):
            # ccs = hapsmash.bamlib.BAM(line)

    chrom2recombination_lst[chrom] = hapsmash_lst 
    chrom2recombination_statistics[chrom] = [
        m.ccs,
        m.lq_ccs,
        m.hq_ccs,
        m.unphaseable_hq_ccs,
        m.phaseable_hq_ccs,
        m.unphased_hq_ccs,
        m.phased_hq_ccs,
        m.hap_consistent_hq_ccs,
        m.hap_inconsistent_hq_ccs,
        m.recombinant_ccs,
    ]
    alignments.close()


def call_recombinantion(
    bam_file: str,
    vcf_file: str,
    region: str,
    region_lst: str,
    min_mapq: int,
    min_alignment_proportion: float,
    min_sequence_identity: float,
    min_bq: int,
    min_trim: float,
    mismatch_window,
    max_mismatch_count: int,
    threads: int,
    version: str,
    out_file: str,
) -> None:

    cpu_start = time.time() / 60
    hapsmash.util.check_caller_input_exists(
        bam_file,
        vcf_file,
        region,
        region_lst,
    )
    _, tname2tsize = hapsmash.bamlib.get_tname2tsize(bam_file)
    chrom_lst, chrom2loci_lst = hapsmash.util.load_loci(region, region_lst, tname2tsize)
    qlen_mean, qlen_lower_limit, qlen_upper_limit, md_threshold = hapsmash.bamlib.get_thresholds(
        bam_file, chrom_lst, tname2tsize
    ) 
    print("hapsmash is calling crossovers and gene conversions with {} threads".format(threads))
    p = mp.Pool(threads)
    manager = mp.Manager()
    chrom2recombination_lst = manager.dict()
    chrom2recombination_statistics = manager.dict()
    get_recombination_arg_lst = [
        (
            chrom,
            tname2tsize[chrom],
            bam_file,
            vcf_file,
            chrom2loci_lst[chrom],
            min_mapq,
            min_alignment_proportion,
            min_sequence_identity,
            qlen_mean,
            qlen_lower_limit, 
            qlen_upper_limit, 
            min_bq,
            min_trim,
            md_threshold,
            mismatch_window,
            max_mismatch_count,
            chrom2recombination_lst,
            chrom2recombination_statistics
        )
        for chrom in chrom_lst
    ]
    p.starmap(
        get_recombination, get_recombination_arg_lst,
    )
    p.close()
    p.join()
    hapsmash.vcflib.dump_recombinantion(chrom_lst, chrom2recombination_lst, out_file)
    hapsmash.vcflib.dump_recombination_statistics(
        chrom_lst, chrom2recombination_statistics, "{}.stat".format(out_file)
    )
    print("hapsmash finished calling crossover and gene conversions")
    cpu_end = time.time() / 60
    duration = cpu_end - cpu_start
    print(
        "hapsmash crossover and gene conversion detection detection took {} minutes".format(
            duration
        )
    )
    hapsmash.util.exit()
