import time
import math
import pysam
import bisect
import natsort
import numpy as np
import himut.bamlib
import hapfusion.util
import hapfusion.bamlib
import hapfusion.haplib
import hapfusion.imglib
import hapfusion.vcflib
import multiprocessing as mp
from dataclasses import dataclass
from collections import defaultdict
from typing import Dict, List, Tuple


@dataclass
class COUNTS:
    ccs: int = 0
    lq_ccs: int = 0 
    hq_ccs: int = 0
    unphased_ccs: int = 0
    hap_consistent_ccs: int = 0 
    hap_inconsistent_ccs: int = 0 
    crossover: int = 0
    gene_conversion: int = 0
    complex_gene_conversion: int = 0


class SMASH:
    def __init__(self):
        self.ccs2mapq = {}
        self.ccs2state = {}
        self.ccs2coord = {}
        self.ccs2hetpos_lst = {} 
        self.ccs2hompos_lst = {}
        self.ccs2denovo_sbs_lst = {}
        self.ccs2denovo_indel_lst = {}
        self.ccs2smash_indel_set = defaultdict(set)
        self.ccs2smash_hetsnp_set = defaultdict(set)
        self.hap2ccs_lst = {hap: [] for hap in ["0", "1", "."]}


def is_low_bq(
    alt: str, 
    min_bq: int, 
    allele2bq_lst: Dict[str, List[int]]
) -> bool:
    alt_bq_lst = allele2bq_lst[himut.util.base2idx[alt]]
    alt_min_bq_count = sum([1 for alt_bq in alt_bq_lst if alt_bq >= min_bq])
    if alt_min_bq_count == 0:
        return True
    else:
        return False


def is_low_gq(
    germ_gq: float,
    min_gq: int,
) -> bool:
    if germ_gq < min_gq:
        return True
    else:
        return False


def is_low_mapq(ccs_mapq: int, min_mapq: int):
    if ccs_mapq < min_mapq:
        return True
    else:
        return False

def is_ccs_phased(hap) -> bool:

    if hap == "0":
        return True
    elif hap == "1":
        return True
    else:
        return False

# def is_chunk(tpos: int, start: int, end: int) -> bool:
#     if start <= tpos and tpos <= end:
#         return True
#     else:
#         return False


# def is_germ_gt(
#     som_gt: str,
#     germ_gt: str,
#     germ_gt_state: str,
#     allelecounts: Dict[int, int],
# ) -> bool:

#     if germ_gt_state == "het":
#         if som_gt == germ_gt:
#             return True
#         else:
#             return False
#     elif germ_gt_state == "hetalt":
#         base_sum = sum(
#             [allelecounts[himut.util.base2idx[base]] for base in himut.util.base_lst]
#         )
#         base_counts = (
#             allelecounts[himut.util.base2idx[germ_gt[0]]]
#             + allelecounts[himut.util.base2idx[germ_gt[1]]]
#         )
#         if base_sum == base_counts and (
#             som_gt[1] == germ_gt[0] or som_gt[1] == germ_gt[1]
#         ):
#             return True
#         else:
#             return False
#     elif germ_gt_state == "homalt":
#         ref_count = allelecounts[himut.util.base2idx[som_gt[0]]]
#         if ref_count == 0 and germ_gt.count(som_gt[1]) == 2:
#             return True
#         else:
#             return False
#     elif germ_gt_state == "homref":
#         if som_gt[1] == germ_gt[0]:
#             return True
#         else:
#             return False


def init_allelecounts():
    rpos2allelecounts = defaultdict(lambda: np.zeros(6))
    rpos2allele2bq_lst = defaultdict(lambda: {0: [], 1: [], 2: [], 3: [], 4: [], 5: []})
    rpos2allele2ccs_lst = defaultdict(
        lambda: {0: [], 1: [], 2: [], 3: [], 4: [], 5: []}
    )
    return rpos2allelecounts, rpos2allele2bq_lst, rpos2allele2ccs_lst, 


def update_allelecounts(
    ccs,
    rpos2allelecounts: Dict[int, np.ndarray],
    rpos2allele2bq_lst: Dict[int, Dict[int, List[int]]],
    rpos2allele2ccs_lst: Dict[int, Dict[int, List[str]]],
):

    tpos = ccs.tstart
    qpos = ccs.qstart
    for (state, ref, alt, ref_len, alt_len) in ccs.cstuple_lst:
        if state == 1:  # match
            for i, alt_base in enumerate(alt):
                epos = tpos + i
                bidx = himut.util.base2idx[alt_base]
                rpos2allelecounts[epos][bidx] += 1
                rpos2allele2ccs_lst[epos][bidx].append(ccs.qname)
                rpos2allele2bq_lst[epos][bidx].append(ccs.bq_int_lst[qpos + i])
        elif state == 2:  # sub
            bidx = himut.util.base2idx[alt]
            rpos2allelecounts[tpos][bidx] += 1
            rpos2allele2ccs_lst[tpos][bidx].append(ccs.qname)
            rpos2allele2bq_lst[tpos][bidx].append(ccs.bq_int_lst[qpos])
        elif state == 3:  # insertion
            rpos2allelecounts[tpos][4] += 1
        elif state == 4:  # deletion
            for j in range(len(ref[1:])):
                rpos2allelecounts[tpos + j][5] += 1
        tpos += ref_len
        qpos += alt_len


def get_hapfusion(
    bam_file: str,
    vcf_file: str,
    ps2hbit_lst: Dict[str, List[str]],
    ps2hpos_lst: Dict[int, List[int]],
    ps2hetsnp_lst: Dict[int, List[Tuple[int, str, str]]], 
    chunkloci_lst: List[Tuple[str, int, int]],
    min_mapq: int,
    min_sequence_identity: float,
    min_alignment_proportion: float,
    qlen_lower_limit: float,
    qlen_upper_limit: float,
    min_bq: int,
    min_trim: float,
    md_threshold: int,
    mismatch_window: int,
    max_mismatch_count: int,
    chrom2hapfusion_lst: Dict[str, List[List[Tuple[str, int, str, str]]]],
    chrom2hapfusion_statistics: Dict[str, List[int]]
) -> List[Tuple[str, int, str, str, int, int, int, float, float]]:

      
    m = COUNTS()
    seen = set()
    counter = 0
    hapfusion_lst = []
    alignments = pysam.AlignmentFile(bam_file, "rb")
    for (chrom, chunk_start, chunk_end) in chunkloci_lst: ## TODO
        phase_set = str(chunk_start)
        hbit_lst = ps2hbit_lst[phase_set] 
        hpos_lst = ps2hpos_lst[phase_set] 
        hetsnp_lst = ps2hetsnp_lst[phase_set] 
        rpos2allelecounts, rpos2allele2bq_lst, rpos2allele2ccs_lst = init_allelecounts()
        for i in alignments.fetch(chrom, chunk_start, chunk_end):
            m.ccs += 1
            ccs = hapfusion.bamlib.BAM(i)
            if not ccs.is_primary: 
                continue
            update_allelecounts(ccs, rpos2allelecounts, rpos2allele2bq_lst, rpos2allele2ccs_lst)
            if is_low_mapq(ccs.mapq, min_mapq):
                continue
            if not (qlen_lower_limit < ccs.qlen and ccs.qlen < qlen_upper_limit):
                continue

            m.hq_ccs += 1
            hapfusion.haplib.get_ccs_hap(ccs, hbit_lst, hpos_lst, hetsnp_lst)  
            if ccs.hap == ".":
                m.unphased_ccs += 1
                continue
            elif is_ccs_phased(ccs.hap):
                m.hap_consistent_ccs += 1
                continue 
            
            if ccs.qname in seen:
                continue
            seen.add(ccs.qname)
            counter += 1
            if counter > 100:
                break
            # def get_phase_switch()(
#     ccs_hap: str,
#     h0_hbit_lst: List[str], 
#     ccs_hbit_lst: List[str],
#     phased_hetsnp_subset_lst: List[Tuple[int, str, str]]
# ) -> Tuple[List[int], str]:

#     haphunter_hetsnp_lst = []
#     haphunter_indel_idx_set = set()
#     haphunter_hetsnp_idx_set = set()
#     if ccs_hap == "0":
#         for i, (j, k) in enumerate(zip(h0_hbit_lst, ccs_hbit_lst)):
#             if k == "-":
#                 haphunter_indel_idx_set.add(i)
#                 continue
#             if j != k:
#                 haphunter_hetsnp_idx_set.add(i)
#                 haphunter_hetsnp_lst.append(phased_hetsnp_subset_lst[i])
#     elif ccs_hap == "1":
#         for i, (j, k) in enumerate(zip(h0_hbit_lst, ccs_hbit_lst)):
#             if k == "-":
#                 haphunter_indel_idx_set.add(i)
#                 continue
#             if j == k:
#                 haphunter_hetsnp_idx_set.add(i)
#                 haphunter_hetsnp_lst.append(phased_hetsnp_subset_lst[i])
#     else:
#         h0_distance = 0
#         h1_distance = 0
#         for (i, j) in zip(h0_hbit_lst, ccs_hbit_lst):
#             if i != j:
#                 h0_distance += 1
#             else:
#                 h1_distance += 1
#         if h0_distance < h1_distance: # ccs haplotype: 0 and "."
#             for i, (j, k) in enumerate(zip(h0_hbit_lst, ccs_hbit_lst)):
#                 if k == "-":
#                     haphunter_indel_idx_set.add(i)
#                     continue
#                 if j != k:
#                     haphunter_hetsnp_idx_set.add(i)
#                     haphunter_hetsnp_lst.append(phased_hetsnp_subset_lst[i])
#         else: # ccs haplotype: 1 
#             for i, (j, k) in enumerate(zip(h0_hbit_lst, ccs_hbit_lst)):
#                 if k == "-":
#                     haphunter_indel_idx_set.add(i)
#                     continue
#                 if j == k:
#                     haphunter_hetsnp_idx_set.add(i)
#                     haphunter_hetsnp_lst.append(phased_hetsnp_subset_lst[i])
#     return haphunter_hetsnp_lst, haphunter_hetsnp_idx_set, haphunter_indel_idx_set 


            # print(ccs.qname) 
            # print(ccs.hap)
            # print(ccs.hbit) 
            # print(ccs.h0_hbit)
            # print(ccs.h1_hbit)
            # print() 
            # ## get_non_crossover()
            ## get_recombination()
            ## phase_switch_count 
            # ccs_hapfusion_hetsnp_candidate_lst = hapfusion.haplib.get_hapfusion_hetsnps(ccs.hap, h0_hbit_lst, ccs_hbit_lst, phased_hetsnp_subset_lst)[0] # search candidate 
            # if len(ccs_hapfusion_hetsnp_candidate_lst) == 0: ## hetsnp is deleted ## is this the result of MMR or sequencing error?
                # continue
               
            # ccs_hapfusion_hetsnp_lst = []
            # trimmed_qstart = math.floor(min_trim * ccs.qlen)
            # trimmed_qend = math.ceil((1 - min_trim) * ccs.qlen)
            # for (tpos, ref, alt) in ccs_hapfusion_hetsnp_candidate_lst:
            #     qpos = ccs.tpos2qpos[tpos]
            #     _, qbq = ccs.tpos2qbase[tpos]
            #     if qbq < min_bq:
            #         continue

            #     if qpos < trimmed_qstart:
            #         continue
            #     if qpos > trimmed_qend:
            #         continue

            #     mismatch_start, mismatch_end = hapfusion.util.get_mismatch_range(tpos, qpos, ccs.qlen, mismatch_window)
            #     jdx = bisect.bisect_left(ccs.mismatch_tpos_lst, mismatch_start)
            #     kdx = bisect.bisect_right(ccs.mismatch_tpos_lst, mismatch_end)
            #     mismatch_count = kdx - jdx 
            #     if mismatch_count > max_mismatch_count:
            #         continue

            #     ins_count = allelecounts[tpos][4]
            #     del_count = allelecounts[tpos][5]
            #     total_count = sum(allelecounts[tpos])
            #     if del_count != 0 or ins_count != 0:
            #         continue
                
            #     if total_count > md_threshold:
            #         continue
            #     ccs_hapfusion_hetsnp_lst.append((tpos, ref, alt))

    #         if len(ccs_hapfusion_hetsnp_lst) != 0:
    #             tend_lst = []
    #             tstart_lst = []
    #             hetsnp_lst = []
    #             homsnp_lst = []
    #             hetindel_lst = []
    #             smash = SMASH()
    #             if len(ccs_hapfusion_hetsnp_lst) % 2 == 0:
    #                 dpos = ccs_hapfusion_hetsnp_lst[int(len(ccs_hapfusion_hetsnp_lst) / 2)][0]
    #                 upos = ccs_hapfusion_hetsnp_lst[int(len(ccs_hapfusion_hetsnp_lst) / 2) - 1][0]
    #                 smash.pos = int((upos + dpos)/2)
    #             else:
    #                 smash.pos = ccs_hapfusion_hetsnp_lst[int(len(ccs_hapfusion_hetsnp_lst) / 2)][0]
    #                 upos = smash.pos   
    #                 dpos = smash.pos + 1
                    
    #             for j in alignments.fetch(ccs.tname, upos, dpos):
    #                 read = hapfusion.bamlib.BAM(j)
    #                 # if read.is_primary:
    #                 if read.is_primary and read.mapq >= 1: 
    #                 # if read.is_primary and read.mapq >= min_mapq: 
    #                     hapfusion.cslib.cs2tpos2qbase(read)
    #                     read.load_mutations(hetset, homset, hetpos_set, hompos_set, phased_hetsnp_set, pos2allele) 
    #                     read_hbit_lst, h0_hbit_lst, h1_hbit_lst, phased_hetsnp_subset_lst = hapfusion.haplib.get_ccs_hbit_lst(
    #                         read,
    #                         hetsnp2bidx,
    #                         hetsnp2hstate,
    #                         phased_hpos_lst,
    #                         phased_hetsnp_lst,
    #                     )
    #                     if len(read_hbit_lst) == 0: # read belongs to more than one haplotype block
    #                         continue 
    #                     tend_lst.append(read.tend)
    #                     tstart_lst.append(read.tstart)
    #                     read.hap = hapfusion.haplib.get_ccs_haplotype(h0_hbit_lst, read_hbit_lst) 
    #                     smash.ccs2mapq[read.qname] = read.mapq
    #                     smash.hap2ccs_lst[read.hap].append(read.qname)
    #                     smash.ccs2coord[read.qname] = (read.tstart, read.tend)
    #                     homsnp_lst.extend(ccs.homsnp_lst)
    #                     hetindel_lst.extend(ccs.hetindel_lst)
    #                     hetsnp_lst.extend(phased_hetsnp_subset_lst)
    #                     phased_hetpos_subset_lst = [phased_hetsnp[0] for phased_hetsnp in phased_hetsnp_subset_lst]
    #                     if read_hbit_lst == h0_hbit_lst: 
    #                         smash.ccs2state[read.qname] = False
    #                         smash.ccs2hetpos_lst[read.qname] = phased_hetpos_subset_lst
    #                         continue
    #                     elif read_hbit_lst == h1_hbit_lst: 
    #                         smash.ccs2state[read.qname] = False
    #                         smash.ccs2hetpos_lst[read.qname] = phased_hetpos_subset_lst
    #                         continue
    #                     seen.add(read.qname)
    #                     smash.ccs2state[read.qname] = True                           
    #                     smash.ccs2hetpos_lst[read.qname] = phased_hetpos_subset_lst
    #                     smash.ccs2denovo_sbs_lst[read.qname] = ccs.denovo_sbs_lst 
    #                     smash.ccs2denovo_indel_lst[read.qname] = ccs.denovo_indel_lst 
    #                     hapfusion_hetsnp_set, hapfusion_indel_set = hapfusion.haplib.get_hapfusion_hetsnps(read.hap, h0_hbit_lst, read_hbit_lst, phased_hetsnp_subset_lst)[1:]
    #                     smash.ccs2smash_indel_set[read.qname] = hapfusion_indel_set
    #                     smash.ccs2smash_hetsnp_set[read.qname] = hapfusion_hetsnp_set 
    #             smash.homsnp_lst = natsort.natsorted(set(homsnp_lst))
    #             smash.hetsnp_lst = natsort.natsorted(set(hetsnp_lst))
    #             smash.hetindel_lst = natsort.natsorted(set(hetindel_lst))
    #             hapfusion_lst.append(
    #                 [
    #                     smash.pos,
    #                     smash.ccs2mapq,
    #                     smash.ccs2state,
    #                     smash.ccs2coord, 
    #                     dict(smash.ccs2smash_indel_set),
    #                     dict(smash.ccs2smash_hetsnp_set),
    #                     smash.ccs2hetpos_lst,
    #                     smash.ccs2denovo_sbs_lst,
    #                     smash.ccs2denovo_indel_lst,
    #                     dict(smash.hap2ccs_lst),
    #                     smash.homsnp_lst,
    #                     smash.hetsnp_lst,
    #                     smash.hetindel_lst,
    #                     [min(tstart_lst), max(tend_lst)],
    #                 ]
    #             ) 
        

    m.lq_ccs = m.ccs - m.hq_ccs
    m.hap_inconsistent_ccs = m.hq_ccs - m.hap_consistent_ccs
    # chrom2hapfusion_lst[chrom] = hapfusion_lst 
    # chrom2hapfusion_statistics[chrom] = [
    #     m.ccs,
    #     m.lq_ccs,
    #     m.hq_ccs,
    #     m.hap_consistent_ccs,
    #     m.hap_inconsistent_ccs,
    #     m.gene_conversion,
    #   m.complex_gene_conversion,
    #   m.crossover,
    # ]
    alignments.close()


def call_recombinations(
    bam_file: str,
    vcf_file: str,
    region: str,
    region_lst: str,
    min_mapq: int,
    min_sequence_identity: float,
    min_alignment_proportion: float,
    min_bq: int,
    min_trim: float,
    mismatch_window,
    max_mismatch_count: int,
    threads: int,
    version: str,
    pdf_dir: str,
    out_file: str,
) -> None:

    # init
    cpu_start = time.time() / 60
    _, tname2tsize = himut.bamlib.get_tname2tsize(bam_file)
    chrom_lst = himut.util.load_loci(region, region_lst, tname2tsize)[0]
    hapfusion.util.check_caller_input_exists(
        bam_file,
        vcf_file,
        out_file,
        chrom_lst,
        tname2tsize
    )
    (
        chrom2ps2hbit_lst,
        chrom2ps2hpos_lst,
        chrom2ps2hetsnp_lst,
        chrom2chunkloci_lst,
    ) = himut.vcflib.load_phased_hetsnps(vcf_file, chrom_lst, tname2tsize)
    qlen_lower_limit, qlen_upper_limit, md_threshold = himut.bamlib.get_thresholds(bam_file, chrom_lst, tname2tsize)

    # start 
    print("hapfusion is calling crossovers and gene conversions with {} threads".format(threads))
    p = mp.Pool(threads)
    manager = mp.Manager()
    chrom2hapfusion_lst = manager.dict()
    chrom2hapfusion_statistics = manager.dict()
    get_hapfusion_arg_lst = [
        (
            bam_file,
            vcf_file,
            chrom2ps2hbit_lst[chrom],
            chrom2ps2hpos_lst[chrom],
            chrom2ps2hetsnp_lst[chrom], 
            chrom2chunkloci_lst[chrom],
            min_mapq,
            min_sequence_identity,
            min_alignment_proportion,
            qlen_lower_limit, 
            qlen_upper_limit, 
            min_bq,
            min_trim,
            md_threshold,
            mismatch_window,
            max_mismatch_count,
            chrom2hapfusion_lst,
            chrom2hapfusion_statistics
        )
        for chrom in chrom_lst
    ]
    p.starmap(
        get_hapfusion, get_hapfusion_arg_lst,
    )
    p.close()
    p.join()
    # end
    sample = hapfusion.vcflib.get_sample(vcf_file)
    # hapfusion.vcflib.dump_hapfusion(chrom_lst, chrom2hapfusion_lst, out_file)
    hapfusion.imglib.dump_hapfusion_pdf(sample, chrom_lst, chrom2hapfusion_lst, pdf_dir)
    # hapfusion.vcflib.dump_hapfusion_statistics(
    #     chrom_lst, chrom2hapfusion_statistics, "{}.stat".format(out_file)
    # )
    print("hapfusion finished calling crossover and gene conversions")
    cpu_end = time.time() / 60
    duration = cpu_end - cpu_start
    print(
        "hapfusion crossover and gene conversion detection detection took {} minutes".format(
            duration
        )
    )
    hapfusion.util.exit()
