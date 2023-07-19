import time
import math
import pysam
import bisect
import natsort
import numpy as np
import himut.bamlib
import hapfusion.util
import hapfusion.gtlib
import hapfusion.bamlib
import hapfusion.haplib
import hapfusion.imglib
import hapfusion.vcflib
import multiprocessing as mp
from datetime import datetime
from dataclasses import dataclass
from typing import Dict, List, Tuple
from collections import Counter, defaultdict

@dataclass
class LOG:
    ccs: int = 0 # number of unique CCS reads 
    lq_ccs: int = 0 
    hq_ccs: int = 0 # primary alignment CCS reads with that mapq > min_mapq and read length between lower and upper bound thresholds 
    unphased_ccs: int = 0 # CCS reads without hetSNPs
    hap_consistent_ccs: int = 0  # CCS reads with exact matching hetSNPs to haplotype block hetSNPs
    hap_inconsistent_ccs: int = 0 # CCS reads with inexact matching hetSNPs to haplotype block hetSNPs 
    hap_denovo_phase_switch_mutation_ccs: int = 0 # CCS reads where de novo mutation induces a phase switching hetSNP
    hap_recombination_candidate_ccs: int = 0 # meiotic recombination candidates
    ambiguous: int = 0
    crossover: int = 0
    gene_conversion: int = 0
    complex_gene_conversion: int = 0


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


def init_allelecounts():
    rpos2basecounts = defaultdict(lambda: np.zeros(6))
    rpos2base2bq_lst = defaultdict(lambda: {0: [], 1: [], 2: [], 3: [], 4: [], 5: []})
    rpos2base2ccs_lst = defaultdict(
        lambda: {0: [], 1: [], 2: [], 3: [], 4: [], 5: []}
    )
    return rpos2basecounts, rpos2base2bq_lst, rpos2base2ccs_lst, 


def update_basecounts(
    ccs,
    rpos2basecounts: Dict[int, np.ndarray],
    rpos2base2bq_lst: Dict[int, Dict[int, List[int]]],
    rpos2base2ccs_lst: Dict[int, Dict[int, List[str]]],
):

    tpos = ccs.tstart
    qpos = ccs.qstart
    for (state, ref, alt, ref_len, alt_len) in ccs.cstuple_lst:
        if state == 1:  # match
            for i, alt_base in enumerate(alt):
                epos = tpos + i
                bidx = himut.util.base2idx[alt_base]
                rpos2basecounts[epos][bidx] += 1
                rpos2base2ccs_lst[epos][bidx].append(ccs.qname)
                rpos2base2bq_lst[epos][bidx].append(ccs.bq_int_lst[qpos + i])
        elif state == 2:  # sub
            bidx = himut.util.base2idx[alt]
            rpos2basecounts[tpos][bidx] += 1
            rpos2base2ccs_lst[tpos][bidx].append(ccs.qname)
            rpos2base2bq_lst[tpos][bidx].append(ccs.bq_int_lst[qpos])
        elif state == 3:  # insertion
            rpos2basecounts[tpos][4] += 1
        elif state == 4:  # deletion
            for j in range(len(ref[1:])):
                rpos2basecounts[tpos + j][5] += 1
        tpos += ref_len
        qpos += alt_len


def measure_hamming_distance(
    s1: str,
    s2: str
):

    hd = 0               
    hidx_lst = []
    for i, (j, k) in enumerate(zip(s1, s2)):
        if k == "-":
            continue
        elif k == "2":
            continue
        elif j != k:
            hd += 1
            hidx_lst.append(i)
    return hd, hidx_lst


def get_hamming_distance(
    hap: str,
    h0_hbit: str, 
    h1_hbit: str,
    ccs_hbit: str,
):

    if hap == "0'": 
        hd, hidx_lst = measure_hamming_distance(h0_hbit, ccs_hbit)
    elif hap == "1'":
        hd, hidx_lst = measure_hamming_distance(h1_hbit, ccs_hbit)
    elif hap == "2": 
        hd, hidx_lst = measure_hamming_distance(h0_hbit, ccs_hbit)
    return hd, hidx_lst


def get_phase_switch_hetsnps(hd, hidx_lst, hetsnp_lst):

    hetsnp_count = len(hetsnp_lst)
    if hd == 1:
        hidx = hidx_lst[0]    
        hetsnp_subset_lst = [hetsnp_lst[hidx]] ## TODO
    else:
        hidx_start = hidx_lst[0]
        hidx_end = hidx_lst[-1] + 1
        upstream_count = hidx_start
        downstream_count = hetsnp_count - hidx_end
        hetsnp_subset_lst = [hetsnp_lst[i] for i in range(hidx_start, hidx_end)]
        phase_switch_tract_length = hetsnp_subset_lst[-1][0] - hetsnp_subset_lst[0][0] 
    return hetsnp_subset_lst, phase_switch_tract_length


def get_phase_switch_ccs_base(ccs):
 
    ccs_bq_lst = []
    ccs_base_lst = []
    for hetsnp in ccs.hetsnp_lst:
        ccs_base, ccs_bq = ccs.tpos2qbase[hetsnp[0]]
        ccs_bq_lst.append(ccs_bq)
        ccs_base_lst.append(ccs_base)
    return ccs_bq_lst, ccs_base_lst 


def is_recombination_candidate(hd: int):

    if hd == 0:
        return False
    else:
        return True    


def get_phase_switch_count(ccs_lst):

    hetsnp2count = Counter([hetsnp for ccs in ccs_lst for hetsnp in ccs.hetsnp_lst])
    return hetsnp2count


def is_denovo(alt_count: int):

    if alt_count > 1:
        return True
    return False


def is_hetsnp(snp_gt, hetsnp_gt):

    if snp_gt == hetsnp_gt:
        return True
    else:
        return False


def is_low_bq(
    bq: int, 
    min_bq: int, 
):

    if bq < min_bq:
        return True
    else:
        return False


def is_low_gq(
    gq: float,
    min_gq: int,
):

    if gq < min_gq:
        return True
    else:
        return False


def is_trimmed(
    ccs,
    qpos: int,
    min_trim: float
):

    trim_qstart, trim_qend = hapfusion.bamlib.get_trim_range(ccs.qlen, min_trim)
    if qpos < trim_qstart:
        return True
    elif qpos > trim_qend:
        return True
    else:
        return False


def is_mismatch_conflict(
    ccs, 
    tpos: int, 
    qpos: int, 
    mismatch_window: int, 
    max_mismatch_count: int
):
    mismatch_tpos_lst = hapfusion.bamlib.get_mismatch_positions(ccs) 
    mismatch_start, mismatch_end = hapfusion.bamlib.get_mismatch_range(
        tpos, qpos, ccs.qlen, mismatch_window
    )
    mismatch_count = (
        bisect.bisect_right(mismatch_tpos_lst, mismatch_end)
        - bisect.bisect_left(mismatch_tpos_lst, mismatch_start)
        - 1
    )
    if mismatch_count > max_mismatch_count:
        return True
    else:
        return False


def get_nco_recomb_length(
    ccs,
    hidx: int
):
    ucount = hidx
    dcount = len(ccs.hetsnp_lst) - hidx -1
    if ucount > 0 and dcount > 0:
        recomb_state = True
        recomb_length = ccs.hetsnp_lst[hidx+1][0] - ccs.hetsnp_lst[hidx-1][0]
    else:
        recomb_state = False
        recomb_length = "." 
    return recomb_state, recomb_length 


def get_co_cnco_recomb_length(
    ccs,
    hidx_lst: List[int]
):

    upos = ccs.hetsnp_lst[hidx_lst[0]][0]
    dpos = ccs.hetsnp_lst[hidx_lst[-1]][0]
    recomb_length = dpos - upos 
    return recomb_length


def get_weighted_hamming_distance(
    s1: str,
    s2: str,
    bq_lst: List[int],
    idx_lst: List[int]
):
    whd = 0
    for idx in idx_lst:
        b1 = s1[idx]
        b2 = s2[idx]
        b2_epsilon = hapfusion.gtlib.get_epsilon(bq_lst[idx]) 
        if b1 != b2:
            if b2 != 2:
                whd += (1 - b2_epsilon) 
            else:
                whd += b2_epsilon/2 
    return whd    


def is_crossover(
    whd, 
    ehd,
    ucount: int,
    dcount: int,
    recomb_length: int
):

    if (whd > ehd) and (ucount == 0 or dcount == 0) and recomb_length > 500:
        return True
    return False 


def is_noncrossover(
    whd: int, 
    ehd: int,
    ucount: int,
    dcount: int,
    recomb_length: int,
):
    if (whd < ehd) and (ucount > 0 and dcount > 0) and recomb_length > 500:
        return True
    return False


def dump_recombinations(
    bam_file: str,
    vcf_file: str,
    region: str, 
    region_list: str,
    tname2tsize: Dict[str, int],
    min_mapq: int,
    qlen_lower_limit: float,
    qlen_upper_limit: float,
    min_bq: int,
    min_gq: int,
    min_trim: float,
    mismatch_window: int,
    max_mismatch_count: int,
    md_threshold: float,
    germline_snp_prior: float,
    threads: int,
    chrom_lst: List[str], 
    chrom2recomb_lst: Dict[str, List[Tuple[str]]],
    version: str,
    out_file: str,
):

    o = open(out_file, "w") # return
    header_lst = [
        "##fileDate={}".format(datetime.now().strftime("%d%m%Y")),
        "##source=hapfusion",
        "##source_version={}".format(version),
        "##content=hapfusion meiotic recombination",
        "##sample={}".format(hapfusion.vcflib.get_sample(vcf_file)),
        '##FILTER=<ID=PASS,Description="All filters passed">',
        '##FILTER=<ID=LowBQ,Description="Base quality score is below the minimum base quality score of {}">'.format(
            min_bq
        ),
        '##FILTER=<ID=LowGQ,Description="Germline genotype quality score is below the minimum genotype quality score of {}">'.format(
            min_gq
        ),
        '##FILTER=<ID=HetAltSite,Description="HetSite is a HetAltSite">',
        '##FILTER=<ID=HomAltSite,Description="HetSite is a HomAltSite">',
        '##FILTER=<ID=HomRefSite,Description="HetSite is a HomRefSite">',
        '##FILTER=<ID=IndelSite,Description="HetSite intersects an indel site">',
        '##FILTER=<ID=HighDepth,Description="Read depth is above the maximum depth threshold of {:.1f}">'.format(
            md_threshold
        ),
        '##FILTER=<ID=Trimmed,Description="Gene conversion is called near the end of a read">',
        '##FILTER=<ID=MismatchCluster,Description="Gene conversion is adjacent to a mismatch within a given mismatch window">',
        '##FILTER=<ID=Indeterminate,Description="CCS read can be either a crossover or a gene conversion">',
    ]
    
    for tname in natsort.natsorted(list(tname2tsize.keys())):
        header_lst.append(
            "##contig=<ID={},length={}>".format(tname, tname2tsize[tname])
        )

    if region is None and region_list is not None:
        region_param = "--region_list {}".format(region_list)
    elif region is not None and region_list is None:
        region_param = "--region {}".format(region)
    elif region is not None and region_list is not None:
        region_param = "--region_list {}".format(region_list)

    param = "{} --min_mapq {} --qlen_lower_limit {} --qlen_upper_limit {} --min_gq {} --min_bq {} --min_trim {} --mismatch_window {} --max_mismatch_count {} --germline_snv_prior {} --threads {} -o {}".format(
        region_param,
        min_mapq,
        qlen_lower_limit,
        qlen_upper_limit,
        min_gq,
        min_bq,
        min_trim,
        mismatch_window,
        max_mismatch_count,
        germline_snp_prior,
        threads,
        out_file,
    )
    
    cmdline = (
        "##hapfusion_command=hapfusion call -i {} --vcf {} {}".format(
            bam_file, 
            vcf_file,
            param
        )
    )
    header_lst.append(cmdline)
    header_lst.append(
        "#{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format("coord", "qname", "phase_set", "state", "event", "haplotype", "hd/whd", "recombination_length", "ccs_hbit", "h0_hbit", "h1_hbit", "hetsnps", "denovo_mutations")
    )
    
    o.write("{}\n".format("\n".join(header_lst)))
    for chrom in chrom_lst:
        for recomb in chrom2recomb_lst[chrom]:
            o.write("{}\n".format("\t".join(recomb)))
    o.close() 


def dump_recombination_candidates(
    chrom_lst: List[str], 
    chrom2recomb_candidate_lst: Dict[str, List[Tuple[str]]]
):

    o = open("hapfusion_candidates.txt", "w")
    o.write("{}\t{}\t{}\t{}\t{}\t{}\n".format("tcoord", "qname", "hap", "ccs_hbit", "h0_hbit", "h1_hbit"))
    for chrom in chrom_lst:
        for recomb_candidate in chrom2recomb_candidate_lst[chrom]:
            o.write("{}\n".format("\t".join(recomb_candidate)))
    o.close() 



def dump_recombination_statistics(
    chrom_lst: List[str], 
    chrom2recomb_statistics: Dict[str, List[int]]
):
    
    row_names = [
            "ccs",
            "lq_ccs",
            "hq_ccs",
            "unphased_ccs",
            "hap_consistent_ccs",
            "hap_inconsistent_ccs",
            "hap_denovo_phase_switch_mutation_ccs",
            "hap_recombination_candidate_ccs",
            "ambiguous",
            "crossover",
            "simple_gene_conversion",
            "complex_gene_conversion",
        ]

    ncol = len(chrom_lst)
    nrow = len(row_names)
    dt = np.zeros((nrow, ncol))
    for i, chrom in enumerate(chrom_lst):
        for j, count in enumerate(chrom2recomb_statistics[chrom]): 
            dt[j][i] = count
    
    o = open("hapfusion.log", "w")
    col_lst = chrom_lst + ["total"] 
    o.write("{:45}{}\n".format("", "\t".join(col_lst)))
    for k in range(nrow):
        rsum =  str(int(np.sum(dt[k])))
        rlst = [str(int(r)) for r in dt[k].tolist()] + [rsum]
        o.write("{:45}{}\n".format(row_names[k], "\t".join(rlst)))
    o.close()         


def get_recombination(
    bam_file: str,
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
    min_gq: int,
    min_trim: float,
    md_threshold: int,
    mismatch_window: int,
    max_mismatch_count: int,
    germline_snp_prior: int,
    chrom2recomb_lst: Dict[str, List[List[Tuple[str, int, str, str]]]],
    chrom2recomb_statistics: Dict[str, List[int]],
    chrom2recomb_candidate_lst: Dict[str, List[Tuple[str]]]
) -> List[Tuple[str, int, str, str, int, int, int, float, float]]:

      
    counter = 0
    m = LOG()
    recomb_lst = []
    co_seen = set()
    ccs_seen = set()
    nco_seen = set()
    recomb_candidate_lst = []
    hapfusion.gtlib.init(germline_snp_prior)
    alignments = pysam.AlignmentFile(bam_file, "rb")
    for (chrom, chunk_start, chunk_end) in chunkloci_lst: ## TODO
        phase_set = str(chunk_start)
        hbit_lst = ps2hbit_lst[phase_set] 
        hpos_lst = ps2hpos_lst[phase_set] 
        hetsnp_lst = ps2hetsnp_lst[phase_set] 
        recombination_candidate_ccs_lst = [] 
        rpos2basecounts, rpos2base2bq_lst, rpos2base2ccs_lst = init_allelecounts()
        for i in alignments.fetch(chrom, chunk_start, chunk_end): ## collect candidates
            ccs = hapfusion.bamlib.BAM(i)
            if not ccs.is_primary: 
                continue
            
            m.ccs += 1
            update_basecounts(ccs, rpos2basecounts, rpos2base2bq_lst, rpos2base2ccs_lst)
            if is_low_mapq(ccs.mapq, min_mapq):
                continue
            if not (qlen_lower_limit < ccs.qlen and ccs.qlen < qlen_upper_limit):
                continue

            m.hq_ccs += 1 
            ccs.get_cs2tpos2qbase() # o: ccs.rpos2qpos, ccs.tpos2qbase, ccs.mismatch_lst
            hapfusion.haplib.get_ccs_hap(ccs, hbit_lst, hpos_lst, hetsnp_lst)  
            if ccs.hap == ".":
                m.unphased_ccs += 1
                continue
            elif is_ccs_phased(ccs.hap):
                m.hap_consistent_ccs += 1
                continue 
            if ccs.qname in ccs_seen:
                continue

            ccs.hd, ccs.hidx_lst = get_hamming_distance(ccs.hap, ccs.h0_hbit, ccs.h1_hbit, ccs.hbit)
            if not is_recombination_candidate(ccs.hd): # hd = 0
                m.hap_denovo_phase_switch_mutation_ccs += 1 
                continue
            recombination_candidate_ccs_lst.append(ccs)
            ccs.get_tcoord()
           
        hetsnp2phase_switch_count = get_phase_switch_count(recombination_candidate_ccs_lst)
        for ccs in recombination_candidate_ccs_lst: ## iterate through candidates
            if ccs.hd == 1:
                hidx = ccs.hidx_lst[0]
                hetsnp = ccs.hetsnp_lst[hidx]
                if hetsnp in nco_seen:
                    continue

                nco_seen.add(hetsnp)
                phase_switch_count = hetsnp2phase_switch_count[hetsnp]  
                if is_denovo(phase_switch_count): # clonal de novo mutations
                    m.hap_denovo_phase_switch_mutation_ccs += phase_switch_count
                    continue     
                              
                rpos = hetsnp[0] - 1
                m.hap_recombination_candidate_ccs += 1
                hetsnp_gt = "{}{}".format(hetsnp[1], hetsnp[2]) 
                recomb_state, recomb_length = get_nco_recomb_length(ccs, hidx)
                snp_gt, snp_gq, snp_state = hapfusion.gtlib.get_germ_gt(hetsnp[1], rpos2base2bq_lst[rpos])
                mut = ",".join(["{}:{}_{}/{}".format(ccs.tname, pos, ref, alt) for (pos, ref, alt) in ccs.hetsnp_lst])
                recomb_candidate_lst.append(
                    [
                        ccs.tcoord,
                        ccs.qname,
                        ccs.hap,
                        ccs.hbit,
                        ccs.h0_hbit,
                        ccs.h1_hbit,
                    ]                    
                )
                
                if not is_hetsnp(snp_gt, hetsnp_gt):
                    recomb_lst.append(
                        [
                            ccs.tcoord,
                            ccs.qname,
                            phase_set,
                            snp_state,
                            "NCO_candidate",
                            ccs.hap,
                            "1",
                            str(recomb_length),
                            ccs.hbit,
                            ccs.h0_hbit, 
                            ccs.h1_hbit,
                            mut,
                            "."
                        ]
                    )
                    continue
                
                bq = ccs.tpos2qbase[hetsnp[0]][1]
                if is_low_bq(bq, min_bq):
                    recomb_lst.append(
                        [
                            ccs.tcoord,
                            ccs.qname,
                            phase_set,
                            "LowBQ",
                            "NCO_candidate",
                            ccs.hap,
                            "1",
                            str(recomb_length),
                            ccs.hbit,
                            ccs.h0_hbit, 
                            ccs.h1_hbit,
                            mut,
                            "."
                        ]
                    )
                    continue
                
                if is_low_gq(snp_gq, min_gq):
                    recomb_lst.append(
                        [
                            ccs.tcoord,
                            ccs.qname,
                            phase_set,
                            "LowGQ",
                            "NCO_candidate",
                            ccs.hap,
                            "1",
                            str(recomb_length),
                            ccs.hbit,
                            ccs.h0_hbit, 
                            ccs.h1_hbit,
                            mut,
                            "."
                        ]
                    )
                    continue
               
                basecounts = rpos2basecounts[rpos]
                read_depth, indel_count = hapfusion.bamlib.get_read_depth(basecounts)
                if indel_count > 0:
                    recomb_lst.append(
                        [
                            ccs.tcoord,
                            ccs.qname,
                            phase_set,
                            "IndelSite",
                            "NCO_candidate",
                            ccs.hap,
                            "1",
                            str(recomb_length),
                            ccs.hbit,
                            ccs.h0_hbit, 
                            ccs.h1_hbit,
                            mut,
                            "."
                        ]
                    )
                    continue
                
                if read_depth > md_threshold:
                    recomb_lst.append(
                        [
                            ccs.tcoord,
                            ccs.qname,
                            phase_set,
                            "HighDepth",
                            "NCO_candidate",
                            ccs.hap,
                            "1",
                            str(recomb_length),
                            ccs.hbit,
                            ccs.h0_hbit, 
                            ccs.h1_hbit,
                            mut,
                            "."
                        ]
                    )
                    continue
                
                qpos = ccs.rpos2qpos[rpos]
                if is_trimmed(ccs, qpos, min_trim):
                    recomb_lst.append(
                        [
                            ccs.tcoord,
                            ccs.qname,
                            phase_set,
                            "Trimmed",
                            "NCO_candidate",
                            ccs.hap,
                            "1",
                            str(recomb_length),
                            ccs.hbit,
                            ccs.h0_hbit, 
                            ccs.h1_hbit,
                            mut,
                            "."
                        ]
                    )
                    continue
                
                if is_mismatch_conflict(ccs, hetsnp[0], qpos, mismatch_window, max_mismatch_count):
                    recomb_lst.append(
                        [
                            ccs.tcoord,
                            ccs.qname,
                            phase_set,
                            "MismatchCluster",
                            "NCO_candidate",
                            ccs.hap,
                            "1",
                            str(recomb_length),
                            ccs.hbit,
                            ccs.h0_hbit, 
                            ccs.h1_hbit,
                            mut
                        ]
                    )
                    continue
                
                if recomb_state:
                    recomb_lst.append(
                        [
                            ccs.tcoord,
                            ccs.qname, 
                            phase_set,
                            "PASS", 
                            "NCO",
                            ccs.hap,
                            "1",
                            str(recomb_length),
                            ccs.hbit, 
                            ccs.h0_hbit, 
                            ccs.h1_hbit, 
                            mut,
                            "."
                        ]
                    )
                    m.gene_conversion += 1
                else:
                    recomb_lst.append(
                        [
                            ccs.tcoord,
                            ccs.qname, 
                            phase_set,
                            "Indeterminate", 
                            "CO_NCO_candidate",
                            ccs.hap,
                            "1",
                            str(recomb_length),
                            ccs.hbit, 
                            ccs.h0_hbit, 
                            ccs.h1_hbit, 
                            mut,
                            "."
                        ]
                    )
                    m.ambiguous += 1
            else:
                hidx_lst = list(range(ccs.hidx_lst[0], ccs.hidx_lst[-1] + 1))
                ehd = len(hidx_lst) - 1
                ucount = hidx_lst[0]
                dcount = len(ccs.hetsnp_lst) - hidx_lst[-1] - 1  
                recomb_length = get_co_cnco_recomb_length(ccs, hidx_lst)
                ccs_hetsnp_bq_lst, ccs_hetsnp_base_lst = get_phase_switch_ccs_base(ccs)
                mut = ",".join(["{}:{}_{}/{}".format(ccs.tname, pos, ref, alt) for (pos, ref, alt) in ccs.hetsnp_lst])
                recomb_candidate_lst.append(
                    [
                        ccs.tcoord,
                        ccs.qname,
                        ccs.hap,
                        ccs.hbit,
                        ccs.h0_hbit,
                        ccs.h1_hbit,
                    ]                    
                )
                m.hap_recombination_candidate_ccs += 1
                if ccs.hap == "0'":
                    whd = get_weighted_hamming_distance(ccs.h0_hbit, ccs.hbit, ccs_hetsnp_bq_lst, hidx_lst)
                elif ccs.hap == "1'":
                    whd = get_weighted_hamming_distance(ccs.h1_hbit, ccs.hbit, ccs_hetsnp_bq_lst, hidx_lst)
                elif ccs.hap == "2":
                    continue
              

                if is_crossover(whd, ehd, ucount, dcount, recomb_length):
                        recomb_lst.append(
                            [
                                ccs.tcoord,
                                ccs.qname, 
                                phase_set,
                                "PASS", 
                                "CO",
                                ccs.hap,
                                str(whd),
                                str(recomb_length),
                                ccs.hbit, 
                                ccs.h0_hbit, 
                                ccs.h1_hbit, 
                                mut,
                                "."
                            ]
                        )
                        m.crossover += 1
                        continue
                if is_noncrossover(whd, ehd, ucount, dcount, recomb_length):
                        recomb_lst.append(
                            [
                                ccs.tcoord,
                                ccs.qname,
                                phase_set, 
                                "PASS", 
                                "CNCO",
                                ccs.hap,
                                str(whd),
                                str(recomb_length),
                                ccs.hbit, 
                                ccs.h0_hbit, 
                                ccs.h1_hbit, 
                                mut,
                                ","
                            ]
                        )
                        m.complex_gene_conversion += 1
                        continue
                recomb_lst.append(
                    [
                        ccs.tcoord,
                        ccs.qname,
                        phase_set, 
                        "Indeterminate", 
                        "CO_NCO_candidate", 
                        ccs.hap,
                        str(whd),
                        str(recomb_length),
                        ccs.hbit, 
                        ccs.h0_hbit, 
                        ccs.h1_hbit, 
                        mut,
                        ".",
                    ]
                )                  
                m.ambiguous += 1 
            ccs_seen.add(ccs.qname)
        
    m.lq_ccs = m.ccs - m.hq_ccs
    m.hap_inconsistent_ccs = m.hq_ccs - (m.unphased_ccs + m.hap_consistent_ccs) 
    chrom2recomb_lst[chrom] = recomb_lst 
    chrom2recomb_statistics[chrom] = [
        m.ccs,
        m.lq_ccs,
        m.hq_ccs,
        m.unphased_ccs,
        m.hap_consistent_ccs,
        m.hap_inconsistent_ccs,
        m.hap_denovo_phase_switch_mutation_ccs,
        m.hap_recombination_candidate_ccs, 
        m.ambiguous,
        m.crossover,
        m.gene_conversion,
        m.complex_gene_conversion,
    ]
    chrom2recomb_candidate_lst[chrom] = recomb_candidate_lst
    alignments.close()

        
def call_recombinations(
    bam_file: str,
    vcf_file: str,
    region: str,
    region_list: str,
    min_mapq: int,
    min_sequence_identity: float,
    min_alignment_proportion: float,
    min_bq: int,
    min_gq: int,
    min_trim: float,
    mismatch_window,
    max_mismatch_count: int,
    germline_snp_prior: float,
    threads: int,
    version: str,
    out_file: str,
) -> None:

    # init
    cpu_start = time.time() / 60
    tname2tsize = himut.bamlib.get_tname2tsize(bam_file)[1]
    chrom_lst = himut.util.load_loci(region, region_list, tname2tsize)[0]
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
    qlen_lower_limit, qlen_upper_limit, md_threshold = hapfusion.bamlib.get_thresholds(bam_file, chrom_lst, tname2tsize)

    # start 
    print("hapfusion is calling crossovers and gene conversions with {} threads".format(threads))
    p = mp.Pool(threads)
    manager = mp.Manager()
    chrom2recomb_lst = manager.dict()
    chrom2recomb_statistics = manager.dict()
    chrom2recomb_candidate_lst = manager.dict()
    get_recomb_arg_lst = [
        (
            bam_file,
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
            min_gq,
            min_trim,
            md_threshold,
            mismatch_window,
            max_mismatch_count,
            germline_snp_prior,
            chrom2recomb_lst,
            chrom2recomb_statistics,
            chrom2recomb_candidate_lst,
        )
        for chrom in chrom_lst
    ]
    p.starmap(
        get_recombination, get_recomb_arg_lst,
    )
    p.close()
    p.join()
    # end
    
    dump_recombinations(
        bam_file,
        vcf_file,
        region,
        region_list,
        tname2tsize,
        min_mapq,
        qlen_lower_limit,
        qlen_upper_limit,
        min_bq,
        min_gq,
        min_trim,
        mismatch_window,
        max_mismatch_count,
        md_threshold,
        germline_snp_prior,
        threads,
        chrom_lst, 
        chrom2recomb_lst, 
        version,
        out_file
    )
    dump_recombination_statistics(chrom_lst, chrom2recomb_statistics)
    dump_recombination_candidates(chrom_lst, chrom2recomb_candidate_lst)
    print("hapfusion finished calling crossover and gene conversions")
    cpu_end = time.time() / 60
    duration = cpu_end - cpu_start
    print(
        "hapfusion crossover and gene conversion detection detection took {} minutes".format(
            duration
        )
    )
    hapfusion.util.exit()
