#!/usr/bin/env python3

import re
import os
import sys
import math
import time
import pysam
import tabix
import cyvcf2
import random
import natsort
import pyfastx
import argparse
import svgwrite
import numpy as np
import multiprocessing as mp
from svglib.svglib import svg2rlg
from collections import defaultdict
from reportlab.graphics import renderPDF
from typing import Set, Dict, List, Tuple


def parse_args(args):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "-i",
        "--vcf",
        type=str,
        required=True,
        help="BAM file to read",
    )    
    parser.add_argument(
        "--region",
        type=str,
        required=False,
        help="target chromosome",
    )
    parser.add_argument(
        "--region_list",
        type=str,
        required=False,
        help="list of target chromosomes separated by new line"
    )
    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        default=1,
        required=False,
        help="number of threads",
    )
    parser.add_argument(
        "-o",
        "--out",
        type=str,
        required=True,
        help="file to return",
    )
    args = args[1:]
    return parser.parse_args(args)


class VCF:
    def __init__(self, line):
        arr = line.strip().split()
        self.chrom = arr[0]
        self.pos = int(arr[1])
        self.id = arr[2]
        self.ref = arr[3]
        self.alt_lst = arr[4].split(",")
        self.qual = float(arr[5]) if arr[5] != "." else arr[5]
        self.is_pass = True if arr[6] == "PASS" else False
        self.info = arr[7]
        self.format_lst = arr[8].split(":")
        self.sample_format_lst = arr[9].split(":")
        hsh = {i:j for i,j in zip(self.format_lst, self.sample_format_lst)}
        if "GT" in hsh: self.sample_gt = hsh["GT"] 
        if "PS" in hsh: self.sample_phase_set = hsh["PS"] 
        if "AD" in hsh: 
            arr = hsh["AD"].split(",")
            self.ref_count = arr[0]
            self.alt_count_arr = arr[1:]

        self.is_snp = False
        self.is_dbs = False
        self.is_indel = False
        if len(self.alt_lst) == 1:  # bi-allelic
            self.is_biallelic = True
            self.alt = self.alt_lst[0]
            if len(self.ref) == 1 and len(self.alt) == 1:  # snp
                self.is_snp = True
            elif len(self.ref) == len(self.alt) == 2:
                self.is_dbs = True
            elif len(self.ref) > len(self.alt): # del
                self.is_indel = True
            elif len(self.ref) < len(self.alt): # ins
                self.is_indel = True
        else:
            self.is_biallelic = False



def get_sample(vcf_file: str):

    if vcf_file.endswith(".vcf"):
        for line in open(vcf_file).readlines():
            if line.startswith("##"):
                if line.startswith("##contig"):
                    continue
            elif line.startswith("#CHROM"):
                sample = line.strip().split()[-1]
                break
    elif vcf_file.endswith(".bgz"):
        v = cyvcf2.VCF(vcf_file)
        sample = v.samples[0]
    return sample


def get_tname2tsize(vcf_file: str):

    tname2tsize = {}
    if vcf_file.endswith(".vcf"):
        for line in open(vcf_file).readlines():
            if line.startswith("##"):
                if line.startswith("##contig"):
                    arr = line.strip().replace("##contig=<ID=", "").split(",")
                    tname = arr[0]
                    tsize = int(arr[1].replace("length=", "").replace(">", ""))
                    tname2tsize[tname] = tsize
            elif line.startswith("#CHROM"):
                break
    elif vcf_file.endswith(".bgz"):
        for line in cyvcf2.VCF(vcf_file).raw_header.split("\n"):
            if line.startswith("##"):
                if line.startswith("##contig"):
                    arr = line.strip().replace("##contig=<ID=", "").split(",")
                    tname = arr[0]
                    tsize = int(arr[1].replace("length=", "").replace(">", ""))
                    tname2tsize[tname] = tsize
            elif line.startswith("#CHROM"):
                break
    return tname2tsize


def load_chrom_lst(
    region: str, 
    region_file: str, 
) -> Tuple[List[str], List[Tuple[str, int, int]]]:

    chrom_lst = []
    if region is None and region_file is not None:
        for line in open(region_file).readlines():
            chrom_lst.append(line.strip().split()[0])
    elif region is not None and region_file is None:
        chrom_lst = [region]
    elif region is not None and region_file is not None:
        for line in open(region_file).readlines():
            chrom_lst.append(line.strip().split()[0])
    else:
        print("--region or --region_list parameter is required")
        sys.exit()
    return chrom_lst


def load_phased_set(
    vcf_file: str,
    chrom_lst: List[str],
    tname2tsize: Dict[str, int]
) -> Dict[str, List[List[Tuple[int, str]]]]:

    chrom2phase_set2hpos_lst = defaultdict()
    for tname in tname2tsize:
        chrom2phase_set2hpos_lst[tname] = defaultdict(list) 
    if vcf_file.endswith(".vcf"):
        for line in open(vcf_file):
            if line.startswith("#"):
                continue
            v = VCF(line)
            if (v.sample_gt == "0|1" or v.sample_gt == "1|0"):
                chrom2phase_set2hpos_lst[v.chrom][v.sample_phase_set].append(v.pos)
    elif vcf_file.endswith(".bgz"):
        tb = tabix.open(vcf_file)
        for chrom in chrom_lst:
            records = tb.query(chrom, 0, tname2tsize[chrom])
            for record in records:
                v = VCF("\t".join(record))
                if v.sample_gt == "0|1" or v.sample_gt == "1|0":
                    chrom2phase_set2hpos_lst[v.chrom][v.sample_phase_set].append(v.pos)

    chrom_set = set(chrom_lst)
    for tname in tname2tsize:
        if tname not in chrom_set:
            del chrom2phase_set2hpos_lst[tname] 
            continue
    return chrom2phase_set2hpos_lst

def get_hapsum(
    sample: str,
    vcf_file: str,
    chrom_lst: List[str],
    tname2tsize: Dict[str, int]
    # chrom2hapsum: Dict[str, int]
): 
    hapsum = 0    
    chrom2phase_set2hpos_lst = load_phased_set(vcf_file, chrom_lst, tname2tsize)
    for chrom in chrom_lst:
        for phase_set in chrom2phase_set2hpos_lst[chrom]:
            hpos_lst = chrom2phase_set2hpos_lst[chrom][phase_set]
            hapsum += (hpos_lst[-1] - hpos_lst[0])
    print(sample, hapsum) 


def dump_hapsum(
    vcf_file: str, 
    region: str,
    region_file: str,
    threads: int,
    out_file: str
): 

    p = mp.Pool(threads)
    manager = mp.Manager()
    chrom2hapsum = manager.dict()
    sample = get_sample(vcf_file)
    tname2tsize = get_tname2tsize(vcf_file)
    chrom_lst = load_chrom_lst(region, region_file)
    get_hapsum(sample, vcf_file, chrom_lst, tname2tsize)
    # get_hapsum_arg_lst = [
    #     (
    #         chrom,
    #         vcf_file,
    #         chrom2hapsum,
    #     )
    #     for chrom in chrom_lst
    # ]
    # p.starmap(
    #     get_hapsum, get_hapsum_arg_lst,
    # )
    # p.close()
    # p.join()

 
def main():
    options = parse_args(sys.argv)
    dump_hapsum(
        options.vcf, 
        options.region, 
        options.region_list,
        options.threads, 
        options.out
    )
    sys.exit(0)


if __name__ == "__main__":
    main()
