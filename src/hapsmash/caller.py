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
import hapsmash.imglib
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


class SMASH:
    def __init__(self):
        self.ccs2mapq = {}
        self.ccs2state = {}
        self.ccs2coord = {}
        self.ccs2hetpos_lst = {} 
        self.ccs2hompos_lst = {}
        self.ccs2denovo_sbs_lst = {}
        self.ccs2denovo_indel_lst = {}
        self.ccs2smash_set = defaultdict(set)
        self.hap2ccs_lst = {hap: [] for hap in ["0", "1", "."]}


def get_hapsmash(
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
    chrom2hapsmash_lst: Dict[str, List[List[Tuple[str, int, str, str]]]],
    chrom2hapsmash_statistics: Dict[str, List[int]]
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
                if ccs.qname in seen:
                    continue
                
                ccs.hap = hapsmash.haplib.get_ccs_haplotype(h0_hbit_lst, ccs_hbit_lst) 
                ccs_hapsmash_hetsnp_candidate_lst, _ = hapsmash.haplib.get_hapsmash_hetsnps(ccs.hap, h0_hbit_lst, ccs_hbit_lst, phased_hetsnp_subset_lst) # search candidate 
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
                    tend_lst = []
                    tstart_lst = []
                    hetsnp_lst = []
                    homsnp_lst = []
                    hetindel_lst = []
                    smash = SMASH()
                    if len(ccs_hapsmash_hetsnp_lst) % 2 == 0:
                        dpos = ccs_hapsmash_hetsnp_lst[int(len(ccs_hapsmash_hetsnp_lst) / 2)][0]
                        upos = ccs_hapsmash_hetsnp_lst[int(len(ccs_hapsmash_hetsnp_lst) / 2) - 1][0]
                        smash.pos = (upos + dpos)/2 
                    else:
                        smash.pos = ccs_hapsmash_hetsnp_lst[int(len(ccs_hapsmash_hetsnp_lst) / 2)][0]

                    for j in alignments.fetch(ccs.tname, ccs.tstart, ccs.tend):
                        read = hapsmash.bamlib.BAM(j)
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
                            if len(read_hbit_lst) == 0: # read belongs to more than one haplotype block
                                continue 
                            tend_lst.append(read.tend)
                            tstart_lst.append(read.tstart)
                            read.hap = hapsmash.haplib.get_ccs_haplotype(h0_hbit_lst, read_hbit_lst) 
                            smash.ccs2mapq[read.qname] = read.mapq
                            smash.hap2ccs_lst[read.hap].append(read.qname)
                            smash.ccs2coord[read.qname] = (read.tstart, read.tend)
                            homsnp_lst.extend(ccs.homsnp_lst)
                            hetindel_lst.extend(ccs.hetindel_lst)
                            hetsnp_lst.extend(phased_hetsnp_subset_lst)
                            phased_hetpos_subset_lst = [phased_hetsnp[0] for phased_hetsnp in phased_hetsnp_subset_lst]
                            if read_hbit_lst == h0_hbit_lst: 
                                smash.ccs2state[read.qname] = False
                                smash.ccs2hetpos_lst[read.qname] = phased_hetpos_subset_lst
                                continue
                            elif read_hbit_lst == h1_hbit_lst: 
                                smash.ccs2state[read.qname] = False
                                smash.ccs2hetpos_lst[read.qname] = phased_hetpos_subset_lst
                                continue 
                            seen.add(read.qname)
                            smash.ccs2state[read.qname] = True                           
                            smash.ccs2hetpos_lst[read.qname] = phased_hetpos_subset_lst
                            smash.ccs2denovo_sbs_lst[read.qname] = ccs.denovo_sbs_lst 
                            smash.ccs2denovo_indel_lst[read.qname] = ccs.denovo_indel_lst 
                            smash.ccs2smash_set[read.qname] = hapsmash.haplib.get_hapsmash_hetsnps(read.hap, h0_hbit_lst, read_hbit_lst, phased_hetsnp_subset_lst)[1]
                    smash.homsnp_lst = natsort.natsorted(set(homsnp_lst))
                    smash.hetsnp_lst = natsort.natsorted(set(hetsnp_lst))
                    smash.hetindel_lst = natsort.natsorted(set(hetindel_lst))
                    hapsmash_lst.append(
                        [
                            smash.pos,
                            smash.ccs2mapq,
                            smash.ccs2state,
                            smash.ccs2coord, 
                            dict(smash.ccs2smash_set),
                            smash.ccs2hetpos_lst,
                            smash.ccs2denovo_sbs_lst,
                            smash.ccs2denovo_indel_lst,
                            dict(smash.hap2ccs_lst),
                            smash.homsnp_lst,
                            smash.hetsnp_lst,
                            smash.hetindel_lst,
                            [min(tstart_lst), max(tend_lst)],
                        ]
                    ) 
        

    chrom2hapsmash_lst[chrom] = hapsmash_lst 
    chrom2hapsmash_statistics[chrom] = [
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
    chrom2hapsmash_lst = manager.dict()
    chrom2hapsmash_statistics = manager.dict()
    get_hapsmash_arg_lst = [
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
            chrom2hapsmash_lst,
            chrom2hapsmash_statistics
        )
        for chrom in chrom_lst
    ]
    p.starmap(
        get_hapsmash, get_hapsmash_arg_lst,
    )
    p.close()
    p.join()
    sample = hapsmash.vcflib.get_sample(vcf_file)
    hapsmash.imglib.dump_hapsmash_pdf(sample, chrom_lst, chrom2hapsmash_lst)
    # hapsmash.vcflib.dump_hapsmash(chrom_lst, chrom2hapsmash_lst, out_file)
    # hapsmash.vcflib.dump_hapsmash_statistics(
    #     chrom_lst, chrom2hapsmash_statistics, "{}.stat".format(out_file)
    # )
    print("hapsmash finished calling crossover and gene conversions")
    cpu_end = time.time() / 60
    duration = cpu_end - cpu_start
    print(
        "hapsmash crossover and gene conversion detection detection took {} minutes".format(
            duration
        )
    )
    hapsmash.util.exit()
