import os
import json
import time
import math
import pysam
import bisect
import natsort
import numpy as np
import hapmash.util
import hapmash.cslib
import hapmash.bamlib
import hapmash.haplib
import hapmash.vcflib
import multiprocessing as mp
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
    recombination_lst = []
    alignments = pysam.AlignmentFile(bam_file, "rb")
    if os.path.exists(vcf_file) and vcf_file.endswith(".vcf"):
        het_set, hom_set = hapmash.vcflib.load_snp_indel(chrom, vcf_file)

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
    ) = hapmash.vcflib.get_phased_hetsnps(vcf_file, chrom, chrom_len)
    # del hidx2bidx
    # del hblock_lst
    # del hidx2hstate
    # del hetsnp2hidx
    # del hidx2hetsnp
   
    counter = 0
    # ccs2hapmash = {}
    hetsnp2cnt = defaultdict(lambda: 0)
    for loci in loci_lst:
        chunkloci_lst = hapmash.util.chunkloci(loci)
        for chunkloci in chunkloci_lst:
            chunk_start, chunk_end = chunkloci[1:]
            allelecounts = defaultdict(lambda: np.zeros(6))
            if vcf_file.endswith(".bgz"):
                het_set, hom_set = hapmash.vcflib.load_bgz_snp_indel((chrom, chunk_start - qlen_mean, chunk_end + qlen_mean), vcf_file)
            
            for line in alignments.fetch(*chunkloci):
                m.ccs += 1
                ccs = hapmash.bamlib.BAM(line)
                hapmash.util.update_allelecounts(ccs, allelecounts)
                
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
                
                if hapmash.util.get_blast_sequence_identity(ccs) < min_sequence_identity:
                    m.lq_ccs += 1
                    continue
                
                m.hq_ccs += 1
                ccs.load_mutations(het_set, hom_set, phased_hetsnp_set) 
                if len(ccs.hetsnp_lst) == 0: # homozygous region 
                    m.unphaseable_hq_ccs += 1 
                    continue 
                elif len(ccs.hetsnp_lst) == 1: 
                    m.phaseable_hq_ccs += 1
                    continue
                else: 
                    m.phaseable_hq_ccs += 1 
                    
                ccs_hbit_lst, h0_hbit_lst, h1_hbit_lst, phased_hetsnp_subset_lst = hapmash.haplib.get_ccs_hbit_lst(
                    ccs,
                    hetsnp2bidx,
                    hetsnp2hstate,
                    phased_hpos_lst,
                    phased_hetsnp_lst,
                )
                if len(ccs_hbit_lst) == 0: # ccs belongs to more than one haplotype block
                    m.unphased_hq_ccs += 1 
                    continue 
                else: # ccs belongs to one haplotype block
                    m.num_phased_hq_ccs += 1 
                    if ccs_hbit_lst == h0_hbit_lst: 
                        m.hap_consistent_hq_ccs += 1    
                        continue
                    elif ccs_hbit_lst == h1_hbit_lst: 
                        m.hap_consistent_hq_ccs += 1    
                        continue 

                ccs_hapswitch_hetsnp_lst = []
                m.hap_inconsistent_hq_ccs += 1
                ccs_hap = hapmash.haplib.get_ccs_haplotype(h0_hbit_lst, ccs_hbit_lst) 
                ccs_hapswitch_hetsnp_candidate_lst = hapmash.haplib.get_hapswitch_hetsnps(ccs_hap, h0_hbit_lst, ccs_hbit_lst, phased_hetsnp_subset_lst) # search candidate 
                if len(ccs_hapswitch_hetsnp_candidate_lst) == 0: ## shouldn't happen?
                    print("WTF", ccs_hapswitch_hetsnp_candidate_lst)
                    continue
                   
                trimmed_qstart = math.floor(min_trim * ccs.qlen)
                trimmed_qend = math.ceil((1 - min_trim) * ccs.qlen)
                for (_, tpos, ref, alt) in ccs_hapswitch_hetsnp_lst:
                    
                    qpos = ccs.tpos2qpos[tpos]
                    _, qbq = ccs.tpos2qbase[tpos]
                    if qbq < min_bq:
                        continue

                    if qpos < trimmed_qstart:
                        continue
                    if qpos > trimmed_qend:
                        continue

                    mismatch_start, mismatch_end = hapmash.util.get_mismatch_range(tpos, qpos, ccs.qlen, mismatch_window)
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

                    hetsnp2cnt[chrom, tpos, ref, alt] += 1
                    ccs_hapswitch_hetsnp_lst.append([chrom, tpos, ref, alt])
                if len(ccs_hapswitch_hetsnp_lst) == 0:
                    continue


            
                # wmark_lst = []
                # wmark_hetsnp_lst = []
                # hapmash_hetsnp_lst = []
                # for jdx, hetsnp in enumerate(read_phased_hetsnp_lst):
                #     hetsnp = "{0[0]}:{0[1]}_{0[2]}/{0[3]}".format(hetsnp)
                #     if jdx in filtered_idx_lst:
                #         wmark_lst.append("*")
                #         if jdx == 0:
                #             downstream_hetsnp = "{0[0]}:{0[1]}_{0[2]}/{0[3]}".format(read_phased_hetsnp_lst[jdx+1])
                #             wmark_hetsnp_lst.append(";".join([hetsnp, downstream_hetsnp]))
                #         elif jdx == read_phased_hetsnp_cnt - 1:
                #             upstream_hetsnp = "{0[0]}:{0[1]}_{0[2]}/{0[3]}".format(read_phased_hetsnp_lst[jdx-1])
                #             wmark_hetsnp_lst.append(";".join([upstream_hetsnp, hetsnp]))
                #         else:
                #             upstream_hetsnp = "{0[0]}:{0[1]}_{0[2]}/{0[3]}".format(read_phased_hetsnp_lst[jdx-1])
                #             downstream_hetsnp = "{0[0]}:{0[1]}_{0[2]}/{0[3]}".format(read_phased_hetsnp_lst[jdx+1])
                #             wmark_hetsnp_lst.append(";".join([upstream_hetsnp, hetsnp, downstream_hetsnp]))
                #         hapmash_hetsnp_lst.append(hetsnp)
                #         hetsnp2cnt[hetsnp] += 1
                #     else:
                #         wmark_lst.append(" ")
                # ccs2hapmash[read.qname] = [read.qname, read_haplotype, hapmash_hetsnp_lst, ",".join(wmark_hetsnp_lst), "".join(wmark_lst), "".join(h0_hbit_lst), "".join(h1_hbit_lst), "".join(read_hbit_lst)]
            # counter += 1
            # if counter > 3:
            #     break
       
    # for ccs in ccs2hapmash:
    #     _, _, hapmash_hetsnp_lst, _, _, _, _, _ = ccs2hapmash[ccs]
    #     hapmash_hetsnp_cnt = len(hapmash_hetsnp_lst)
    #     if sum([hetsnp2cnt[hetsnp] for hetsnp in hapmash_hetsnp_lst]) == hapmash_hetsnp_cnt:
    #         hapmash_lst.append(ccs2hapmash[ccs])
    #         m.num_hap_chimeric_hq_reads += 1

    chrom2recombination_lst[chrom] = recombination_lst 
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
    hapmash.util.check_caller_input_exists(
        bam_file,
        vcf_file,
        region,
        region_lst,
    )
    _, tname2tsize = hapmash.bamlib.get_tname2tsize(bam_file)
    chrom_lst, chrom2loci_lst = hapmash.util.load_loci(region, region_lst, tname2tsize)
    qlen_mean, qlen_lower_limit, qlen_upper_limit, md_threshold = hapmash.bamlib.get_thresholds(
        bam_file, chrom_lst, tname2tsize
    ) 
    print("hapmash is calling crossovers and gene conversions with {} threads".format(threads))
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
    hapmash.vcflib.dump_recombinantion(chrom_lst, chrom2recombination_lst, out_file)
    hapmash.vcflib.dump_recombination_statistics(
        chrom_lst, chrom2recombination_statistics, "{}.stat".format(out_file)
    )
    print("hapmash finished calling crossover and gene conversions")
    cpu_end = time.time() / 60
    duration = cpu_end - cpu_start
    print(
        "hapmash crossover and gene conversion detection detection took {} minutes".format(
            duration
        )
    )
    hapmash.util.exit()
