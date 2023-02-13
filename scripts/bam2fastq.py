#!/usr/bin/env python3

import sys
import pysam
import argparse
import himut.bamlib
from typing import Dict, List, Tuple, Set


def parse_args(args):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "-i",
        "--bam",
        type=str,
        required=True,
        help="BAM filed to read",
    )
    parser.add_argument(
        "--bed",
        type=str,
        required=True,
        help="list of coordinates separated by new line"
    ) 
    parser.add_argument(
        "--min_mapq",
        type=int,
        default=60,
        required=False,
        help="minimum mapping quality score",
    )
    parser.add_argument(
        "-o",
        "--out",
        type=str,
        required=True,
        help="file to write",
    )
    args = args[1:]
    return parser.parse_args(args)


def load_loci_lst(
    region_file: str, 
) -> Tuple[List[str], List[Tuple[str, int, int]]]:

    loci_lst = []
    for line in open(region_file).readlines():
        chrom, start, end = line.strip().split()
        loci_lst.append((chrom, int(start), int(end)))
    return loci_lst


def dump_ccs_statistics(
    bam_file: str,
    bed_file: str,
    min_mapq: int,
    out_file: str
): 

    o = open(out_file, "w")    
    loci_lst = load_loci_lst(bed_file)
    alignments = pysam.AlignmentFile(bam_file, "rb")
    for loci in loci_lst:
        for line in alignments.fetch(*loci):
            ccs = himut.bamlib.BAM(line)
            if not ccs.is_primary:
                continue
            if ccs.mapq < min_mapq:
                continue
            qbq = "".join([chr(bq+33) for bq in ccs.bq_int_lst])
            o.write("@{}\n{}\n+\n{}\n".format(ccs.qname, ccs.qseq, qbq)) 
    o.close() 


def main():
    options = parse_args(sys.argv)
    dump_ccs_statistics(
        options.bam, 
        options.bed,
        options.min_mapq,
        options.out
    )
    sys.exit(0)


if __name__ == "__main__":
    main()

