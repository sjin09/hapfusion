import os
import sys
import psutil
import natsort
import numpy as np
from collections import defaultdict
from typing import Dict, List, Tuple


base_lst = list("ATGC+-")
base2idx = {base: idx for idx, base in enumerate(base_lst)}
idx2base = {idx: base for idx, base in enumerate(base_lst)}


class NestedDefaultDict(defaultdict):
    def __init__(self, *args, **kwargs):
        super(NestedDefaultDict, self).__init__(NestedDefaultDict, *args, **kwargs)

    def __repr__(self):
        return repr(dict(self))


def exit():
    print("exiting hapsmash")
    sys.exit(0)


def load_loci(
    region: str, 
    region_list: str, 
    tname2tsize: Dict[str, int]
) -> Tuple[List[str], List[Tuple[str, int, int]]]:

    chrom2loci_lst = defaultdict(list)
    if region is None and region_list is not None:
        for line in open(region_list).readlines():
            arr = line.strip().split()
            chrom = arr[0]
            if len(arr) == 1:
                chrom2loci_lst[chrom].append((chrom, 0, tname2tsize[chrom]))
            else:
                chrom2loci_lst[chrom].append((chrom, int(arr[1]), int(arr[2])))
    elif region is not None and region_list is None:
        chrom_lst = [region]
        if region in tname2tsize:
            chrom2loci_lst[region].append((region, 0, tname2tsize[region]))
        else:
            print("{} does not exist in the BAM file".format(region))
            exit()
    elif region is not None and region_list is not None:
        for line in open(region_list).readlines():
            arr = line.strip().split()
            chrom = arr[0]
            if len(arr) == 1:
                chrom2loci_lst[chrom].append((chrom, 0, tname2tsize[chrom]))
            else:
                chrom2loci_lst[chrom].append((chrom, int(arr[1]), int(arr[2])))
    else:
        for tname, tsize in tname2tsize.items():
            chrom2loci_lst[tname].append((tname, 0, tsize))

    chrom_lst = natsort.natsorted(list(chrom2loci_lst.keys()))
    for chrom in chrom_lst:
        chrom2loci_lst[chrom] = natsort.natsorted(chrom2loci_lst[chrom])
    return chrom_lst, chrom2loci_lst


def chunkloci(loci: Tuple[str, int, int]) -> List[Tuple[str, int, int]]:
    chunk_loci_lst = []
    chrom, start, end = loci
    chunk_start_lst = list(range(start, end, 200000))
    for i, chunk_start in enumerate(chunk_start_lst[:-1]):
        chunk_loci_lst.append((chrom, chunk_start, chunk_start_lst[i+1]))
    if (chrom, chunk_start_lst[-1], end) not in chunk_loci_lst:
        chunk_loci_lst.append((chrom, chunk_start_lst[-1], end))
    return chunk_loci_lst


def check_num_threads(thread_count: int):
    system_thread_count = psutil.cpu_count()
    if thread_count > system_thread_count:
        print("System does not have {} number of threads".format(thread_count))
        print("Please provide a more appropriate number of threads")
        exit()


def check_bam_file(bam_file: str):

    if bam_file is None:
        print("Please provide the path to the BAM file")
        return 1
    else:
        if bam_file.endswith(".bam"):
            idxfile = "{}.bai".format(bam_file)
            if os.path.exists(bam_file) and os.path.exists(idxfile):
                if os.path.getsize(bam_file) != 0 and os.path.getsize(idxfile) != 0:
                    return 0
                else:
                    return 1
            elif not os.path.exists(bam_file) and os.path.exists(idxfile):
                print("BAM file is missing")
                return 1
            elif os.path.exists(bam_file) and not os.path.exists(idxfile):
                print("BAM index file is missing")
                print("Use samtools index to index your BAM file")
                return 1
            elif not os.path.exists(bam_file) and not os.path.exists(idxfile):
                print("BAM file is missing")
                print("BAM index file is missing")
                return 1
        else:
           print("Did you provide a BAM file?") 
           print("BAM files have to a .bam suffix") 


def check_vcf_file(vcf_file: str):
    if vcf_file is None:
        print("Please provide the path to the VCF file")
        return 1
    else: 
        if os.path.exists(vcf_file):
            if vcf_file.endswith(".vcf"):
                if os.path.getsize(vcf_file) == 0:
                    return 1
                else:
                    return 0
            elif vcf_file.endswith(".bgz"):
                tbi_file = vcf_file + ".tbi"  
                if os.path.exists(tbi_file):
                    return 0 
                else:
                    print("tabix index file does not exist")
                    return 1
            elif vcf_file.endswith(".gz"):
                print("hapsmash doesn't support loading of gzip compressed VCF files") 
                return 1
            else:
                print("VCF file must have the .vcf suffix")
                return 1
        else:
            print("{} file is missing".format(vcf_file))
            return 1



def check_vcf_content(
    chrom_lst: List[str], 
    chrom2sub: Dict[str, Tuple[str, int, str, str]],
    vcf_file: str,
):

    state = 0
    for chrom in chrom_lst:
        num_var = len(chrom2sub[chrom])
        if num_var == 0:
            state = 1
        print("{}: loading {} number of variants from {}".format(vcf_file, num_var, chrom))
    if state == 1:
        print("VCF file might be corrupted")
        print("Please check the contents of the input VCF files")
        exit()
       

def check_out_file(out_file: str):

    if out_file.endswith(".vcf"):
        return 0 
    elif out_file.endswith(".gz"):
        print("hapsmash doesn't support return gzip compressed VCF files") 
        return 1
    elif out_file.endswith(".bgz"):
        print("hapsmash doesn't support return bgzip compressed VCF files") 
        return 1
    else:
        print("hapsmash doesn't recgonise the suffix of the file")
        return 1

    
def check_caller_input_exists(
    bam_file: str,
    vcf_file: str,
    region: str,
    region_lst: str,
) -> None:

    counter = 0
    counter += check_bam_file(bam_file)
    counter += check_vcf_file(vcf_file)
    if region is None and region_lst is None:
        counter += 1
        print("--region and --region_list params are missing")
        print("Please provide a chromosome or a file with a list of chromosomes")
        
    if counter > 0:
        print("One or more inputs and parameters are missing")
        print("Please provide the correct inputs and parameters")
        exit()


def check_phaser_input_exists(
    bam_file: str,
    vcf_file: str,
    out_file: str, 
) -> None:

    counter = 0
    counter += check_bam_file(bam_file)
    counter += check_vcf_file(vcf_file)
    counter += check_out_file(out_file)
    if not out_file.endswith(".vcf"):
        print("Please use the suffix .vcf for the output file")
        counter += 1

    if counter > 0:
        print("One or more inputs and parameters are missing")
        print("Please provide the correct inputs and parameters")
        exit()


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
    ccs,
    allelecounts: Dict[int, np.ndarray],
):

    tpos = ccs.tstart
    qpos = ccs.qstart
    ccs.tpos2qpos = {}
    ccs.tpos2qbase = {} 
    if ccs.is_primary:
        for cstuple in ccs.cstuple_lst:
            state, ref, alt, ref_len, alt_len, = cstuple
            if state == 1:  # match
                for i, alt_base in enumerate(alt):
                    ccs.tpos2qpos[tpos + i + 1] = qpos + i
                    allelecounts[tpos + i + 1][base2idx[alt_base]] += 1
                    ccs.tpos2qbase[tpos + i + 1] = (alt_base, ccs.bq_int_lst[qpos + i])
            elif state == 2:  # sub
                ccs.tpos2qpos[tpos + 1] = qpos
                allelecounts[tpos + 1][base2idx[alt]] += 1
                ccs.tpos2qbase[tpos + 1] = (alt, ccs.bq_int_lst[qpos])
            elif state == 3:  # insertion
                allelecounts[tpos + 1][4] += 1
            elif state == 4:  # deletion
                for j in range(len(ref[1:])):
                    allelecounts[tpos + j + 1][5] += 1
                    ccs.tpos2qbase[tpos + j + 1] = ("-", 0)
            tpos += ref_len
            qpos += alt_len
    # else:
    #     return hapsmash.cslib.cs2tpos2qbase(ccs)

