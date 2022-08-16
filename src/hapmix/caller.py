import json
import time
import math
import pysam
import bisect
import natsort
import numpy as np
from regex import W
import hapmix.util
import hapmix.cslib
import hapmix.bamlib
import hapmix.haplib
import hapmix.vcflib
import multiprocessing as mp
from dataclasses import dataclass, field
from collections import defaultdict, Counter
from typing import Set, Dict, List, Tuple


@dataclass
class METRICS:
    num_reads: int = 0
    num_lq_reads: int = 0 
    num_hq_reads: int = 0
    num_unphaseable_hq_reads: int = 0
    num_phaseable_hq_reads: int = 0
    num_unphased_hq_reads: int = 0 
    num_phased_hq_reads: int = 0 
    num_hap_consistent_hq_reads: int = 0 
    num_hap_inconsistent_hq_reads: int = 0 
    num_hap_chimeric_hq_reads: int = 0 

def get_blast_sequence_identity(read) -> float:

    mismatch_count = 0
    for cstuple in read.cstuple_lst:
        mstate, _ref, _alt, ref_len, alt_len = cstuple
        if mstate == 1:  # match
            continue
        elif mstate == 2:  # mismatch: snp
            mismatch_count += 1
        elif mstate == 3:  # mismatch: insertion
            mismatch_count += alt_len
        elif mstate == 4:  # mismatch: deletion
            mismatch_count += ref_len
    blast_sequence_identity = (read.target_alignment_length - mismatch_count)/float(read.target_alignment_length)
    return blast_sequence_identity


def get_mismatch_range(
    tpos: int, 
    qpos: int,
    qlen: int, 
    window: int
):
    qstart, qend = [qpos - window, qpos + window]
    if qstart < 0:
        urange = window + qstart
        drange = window + abs(qstart)
    elif qend > qlen:
        urange = window + abs(qend - qlen)
        drange = qlen - qpos
    else:
        urange = window
        drange = window
    tstart = tpos - urange
    tend = tpos + drange
    return tstart, tend


def update_allelecounts(
    read,
    allelecounts: Dict[int, np.ndarray],
):

    tpos = read.tstart
    qpos = read.qstart
    read.tpos2qpos = {}
    read.tpos2qbase = {} 
    hapmix.cslib.cs2tuple(read)
    if read.mapq >= 1 and read.is_primary:
        for cstuple in read.cstuple_lst:
            state, ref, alt, ref_len, alt_len, = cstuple
            if state == 1:  # match
                for i, alt_base in enumerate(alt):
                    read.tpos2qpos[tpos + i + 1] = qpos + i
                    allelecounts[tpos + i + 1][hapmix.util.base2idx[alt_base]] += 1
                    read.tpos2qbase[tpos + i + 1] = (alt_base, read.bq_int_lst[qpos + i])
            elif state == 2:  # sub
                read.tpos2qpos[tpos + 1] = qpos
                allelecounts[tpos + 1][hapmix.util.base2idx[alt]] += 1
                read.tpos2qbase[tpos + 1] = (alt, read.bq_int_lst[qpos])
            elif state == 3:  # insertion
                allelecounts[tpos + 1][4] += 1
            elif state == 4:  # deletion
                for j in range(len(ref[1:])):
                    allelecounts[tpos + j + 1][5] += 1
                    read.tpos2qbase[tpos + j + 1] = ("-", 0)
            tpos += ref_len
            qpos += alt_len
    else:
        hapmix.cslib.cs2tpos2qbase(read)


def get_hapmix_molecules(
    chrom: str,
    chrom_len: int,
    bam_file: str,
    vcf_file: str,
    loci_lst: List[Tuple[str, int, int]],
    min_mapq: int,
    min_alignment_proportion: float,
    min_sequence_identity: float,
    qlen_lower_limit: float,
    qlen_upper_limit: float,
    min_bq: int,
    min_trim: float,
    md_threshold: int,
    mismatch_window: int,
    max_mismatch_count: int,
    chrom2hapmix_lst: Dict[str, List[List[Tuple[str, int, str, str]]]],
    chrom2hetsnp_statistics: Dict[str, List[int]]
) -> List[Tuple[str, int, str, str, int, int, int, float, float]]:

    # counter = 0
    m = METRICS()
    hapmix_lst = []
    alignments = pysam.AlignmentFile(bam_file, "rb")
    (
        snp_set,
        hpos_lst,
        hetsnp_lst,
        hblock_lst,
        hidx2bidx,
        hidx2hetsnp,
        hidx2hstate,
        hetsnp2bidx,
        hetsnp2hidx,
        hetsnp2hstate,
    ) = hapmix.vcflib.get_germline_mutations(vcf_file, chrom, chrom_len)
    del hblock_lst
    del hidx2hetsnp
    
    counter = 0
    ccs2hapmix = {}
    hetsnp2cnt = defaultdict(lambda: 0)
    for loci in loci_lst:
        chunkloci_lst = hapmix.util.chunkloci(loci)
        for chunkloci in chunkloci_lst:
            allelecounts = defaultdict(lambda: np.zeros(6))
            for line in alignments.fetch(*chunkloci):
                m.num_reads += 1
                read = hapmix.bamlib.BAM(line)
                update_allelecounts(read, allelecounts)
                
                if not read.is_primary:
                    m.num_lq_reads += 1
                    continue

                if read.mapq < min_mapq:
                    m.num_lq_reads += 1
                    continue

                if read.qlen < qlen_lower_limit or read.qlen > qlen_upper_limit:
                    m.num_lq_reads += 1
                    continue
               
                if read.query_alignment_proportion < min_alignment_proportion:
                    m.num_lq_reads += 1
                    continue
                
                if get_blast_sequence_identity(read) < min_sequence_identity:
                    m.num_lq_reads += 1
                    continue
                
                m.num_hq_reads += 1
                read_phased_hetsnp_cnt, read_phased_hetsnp_lst = hapmix.haplib.get_read_hetsnps(read.tpos2qbase, hpos_lst, hetsnp_lst) 
                if read_phased_hetsnp_cnt == 0:
                    m.num_unphaseable_hq_reads += 1 
                    continue 
                else:
                    m.num_phaseable_hq_reads += 1 
                    
                read_hbit_lst, h0_hbit_lst, h1_hbit_lst = hapmix.haplib.get_read_hbit_lst(
                    read.tpos2qbase,
                    hetsnp2bidx,
                    hetsnp2hstate,
                    read_phased_hetsnp_lst,
                )
                if len(read_hbit_lst) == 0:
                    m.num_unphased_hq_reads += 1 
                    continue 
                else:
                    m.num_phased_hq_reads += 1 
                
                if read_hbit_lst == h0_hbit_lst:
                    m.num_hap_consistent_hq_reads += 1    
                    continue
                elif read_hbit_lst == h1_hbit_lst:
                    m.num_hap_consistent_hq_reads += 1    
                    continue 

                idx_lst, read_haplotype = hapmix.haplib.get_hapmix_idx(h0_hbit_lst, read_hbit_lst)
                if len(idx_lst) == 0:
                    continue
                   
                filtered_idx_lst = []
                hapmix.cslib.cs2mismatch(read)
                m.num_hap_inconsistent_hq_reads += 1
                trimmed_qstart = math.floor(min_trim * read.qlen)
                trimmed_qend = math.ceil((1 - min_trim) * read.qlen)
                mismatch_lst = natsort.natsorted([mismatch for mismatch in read.mismatch_lst if mismatch not in snp_set])
                mpos_lst = [mismatch[1] for mismatch in mismatch_lst]
                for idx in idx_lst:
                    tpos = read_phased_hetsnp_lst[idx][1]
                    qbase, qbq = read.tpos2qbase[tpos]
                    if qbase == "-":
                        continue

                    if qbq < min_bq:
                        continue
                    
                    qpos = read.tpos2qpos[tpos]
                    if qpos < trimmed_qstart:
                        continue
                    if qpos > trimmed_qend:
                        continue

                    mismatch_start, mismatch_end = get_mismatch_range(tpos, qpos, read.qlen, mismatch_window)
                    jdx = bisect.bisect_left(mpos_lst, mismatch_start)
                    kdx = bisect.bisect_right(mpos_lst, mismatch_end)
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

                    filtered_idx_lst.append(idx)
                if len(filtered_idx_lst) == 0:
                    continue
            
                wmark_lst = []
                wmark_hetsnp_lst = []
                hapmix_hetsnp_lst = []
                for jdx, hetsnp in enumerate(read_phased_hetsnp_lst):
                    hetsnp = "{0[0]}:{0[1]}_{0[2]}/{0[3]}".format(hetsnp)
                    if jdx in filtered_idx_lst:
                        wmark_lst.append("*")
                        if jdx == 0:
                            downstream_hetsnp = "{0[0]}:{0[1]}_{0[2]}/{0[3]}".format(read_phased_hetsnp_lst[jdx+1])
                            wmark_hetsnp_lst.append(";".join([hetsnp, downstream_hetsnp]))
                        elif jdx == read_phased_hetsnp_cnt - 1:
                            upstream_hetsnp = "{0[0]}:{0[1]}_{0[2]}/{0[3]}".format(read_phased_hetsnp_lst[jdx-1])
                            wmark_hetsnp_lst.append(";".join([upstream_hetsnp, hetsnp]))
                        else:
                            upstream_hetsnp = "{0[0]}:{0[1]}_{0[2]}/{0[3]}".format(read_phased_hetsnp_lst[jdx-1])
                            downstream_hetsnp = "{0[0]}:{0[1]}_{0[2]}/{0[3]}".format(read_phased_hetsnp_lst[jdx+1])
                            wmark_hetsnp_lst.append(";".join([upstream_hetsnp, hetsnp, downstream_hetsnp]))
                        hapmix_hetsnp_lst.append(hetsnp)
                        hetsnp2cnt[hetsnp] += 1
                    else:
                        wmark_lst.append(" ")
                ccs2hapmix[read.qname] = [read.qname, read_haplotype, hapmix_hetsnp_lst, ",".join(wmark_hetsnp_lst), "".join(wmark_lst), "".join(h0_hbit_lst), "".join(h1_hbit_lst), "".join(read_hbit_lst)]
            # counter += 1
            # if counter > 3:
            #     break
       
    for ccs in ccs2hapmix:
        _, _, hapmix_hetsnp_lst, _, _, _, _, _ = ccs2hapmix[ccs]
        hapmix_hetsnp_cnt = len(hapmix_hetsnp_lst)
        if sum([hetsnp2cnt[hetsnp] for hetsnp in hapmix_hetsnp_lst]) == hapmix_hetsnp_cnt:
            hapmix_lst.append(ccs2hapmix[ccs])
            m.num_hap_chimeric_hq_reads += 1

    chrom2hapmix_lst[chrom] = hapmix_lst 
    chrom2hetsnp_statistics[chrom] = [
        m.num_reads,
        m.num_lq_reads,
        m.num_hq_reads,
        m.num_unphaseable_hq_reads,
        m.num_phaseable_hq_reads,
        m.num_unphased_hq_reads,
        m.num_phased_hq_reads,
        m.num_hap_consistent_hq_reads,
        m.num_hap_inconsistent_hq_reads,
        m.num_hap_chimeric_hq_reads,
    ]
    alignments.close()


def call_hapmix(
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
    hapmix.util.check_caller_input_exists(
        bam_file,
        vcf_file,
        region,
        region_lst,
    )
    _, tname2tsize = hapmix.bamlib.get_tname2tsize(bam_file)
    chrom_lst, chrom2loci_lst = hapmix.util.load_loci(region, region_lst, tname2tsize)
    qlen_lower_limit, qlen_upper_limit, md_threshold = hapmix.bamlib.get_thresholds(
        bam_file, chrom_lst, tname2tsize
    ) 
    print("hapmix is calling crossovers and gene conversions with {} threads".format(threads))
    p = mp.Pool(threads)
    manager = mp.Manager()
    chrom2hapmix_lst = manager.dict()
    chrom2hapmix_statistics = manager.dict()
    get_chrom_hapmix_arg_lst = [
        (
            chrom,
            tname2tsize[chrom],
            bam_file,
            vcf_file,
            chrom2loci_lst[chrom],
            min_mapq,
            min_alignment_proportion,
            min_sequence_identity,
            qlen_lower_limit, 
            qlen_upper_limit, 
            min_bq,
            min_trim,
            md_threshold,
            mismatch_window,
            max_mismatch_count,
            chrom2hapmix_lst,
            chrom2hapmix_statistics
        )
        for chrom in chrom_lst
    ]
    p.starmap(
        get_hapmix_molecules, get_chrom_hapmix_arg_lst,
    )
    p.close()
    p.join()
    hapmix.vcflib.dump_recombinantion(chrom_lst, chrom2hapmix_lst, out_file)
    hapmix.vcflib.dump_hapmix_statistics(
        chrom_lst, chrom2hapmix_statistics, "{}.stat".format(out_file)
    )
    print("hapmix finished calling crossover and gene conversions")
    cpu_end = time.time() / 60
    duration = cpu_end - cpu_start
    print(
        "hapmix crossover and gene conversion detection detection took {} minutes".format(
            duration
        )
    )
    hapmix.util.exit()