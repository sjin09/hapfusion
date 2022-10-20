#!/usr/bin/env python3

import re
import sys
import pysam
import tabix
import cyvcf2
import natsort
import argparse
import multiprocessing as mp
from collections import defaultdict
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
        "--region_list",
        type=str,
        required=True,
        help="list of loci separated by new line"
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


class BAM:
    def __init__(self, line):
        # target
        self.tname = line.reference_name
        self.tstart = line.reference_start
        self.tend = line.reference_end
        self.tcoord = "{}:{}-{}".format(self.tname, self.tstart, self.tend)
        # query
        self.qname = line.query_name.replace("/ccs", "")
        self.qstart = line.query_alignment_start
        self.qend = line.query_alignment_end
        self.qseq = line.query_sequence
        # self.qlen = len(self.qseq)
        self.mapq = line.mapping_quality
        self.bq_int_lst = line.query_qualities
        # self.qv = sum(self.bq_int_lst)/self.qlen
        # self.hbq_proportion = self.bq_int_lst.count(93)/float(self.qlen)
        self.query_alignment_length = self.qend - self.qstart
        # self.query_alignment_proportion = self.query_alignment_length/float(self.qlen)
        self.cs_tag = line.get_tag("cs") if line.has_tag("cs") else "."
        if line.has_tag("tp"):
            if line.get_tag("tp") == "P":
                self.is_primary = True
            else:
                self.is_primary = False
        else:
            self.is_primary = False


    def cs2lst(self):
        cslst = [cs for cs in re.split("(:[0-9]+|\*[a-z][a-z]|[=\+\-][A-Za-z]+)", self.cs_tag)]
        cslst = [cs.upper() for cs in cslst if cs != ""]
        return cslst

            
    def cs2tuple(self) -> List[Tuple[int, str, str, int, int]]:
        qpos = self.qstart
        self.cstuple_lst = []
        cs_lst = self.cs2lst()
        for cs in cs_lst:
            m = cs[1:]
            mlen = len(m)
            qstart = qpos
            if cs.startswith("="):  # match # --cs=long
                cs = ":{}".format(mlen)
                t = (1, m, m, mlen, mlen)
            elif cs.startswith(":"):  # match # --cs=short
                mlen = int(m)
                qend = qpos + mlen
                m = self.qseq[qstart:qend]
                t = (1, m, m, mlen, mlen)
            elif cs.startswith("*"):  # snp # target and query
                mlen = 1
                ref, alt = list(m)
                t = (2, ref, alt, 1, 1)
            elif cs.startswith("+"):  # insertion # query
                ref = self.qseq[qpos - 1]
                alt = ref + m
                t = (3, ref, alt, 0, mlen)
            elif cs.startswith("-"):  # deletion # target
                alt = self.qseq[qpos - 1]
                ref = alt + m
                t = (4, ref, alt, mlen, 0)
                mlen = 0
            self.cstuple_lst.append(t)
            qpos += mlen


    def cs2subindel(self):
        self.cs2tuple()
        tpos = self.tstart
        qpos = self.qstart
        self.tsbs_lst = []
        self.qsbs_lst = []
        self.qsbs_bq_lst = []
        self.tindel_lst = []
        for cstuple in self.cstuple_lst:
            state, ref, alt, ref_len, alt_len, = cstuple
            if state == 2:  # snp 
                if ref == "N": 
                    continue
                self.qsbs_bq_lst.append(self.bq_int_lst[qpos])
                self.tsbs_lst.append((tpos + 1, ref, alt))
                self.qsbs_lst.append((qpos + 1, alt, ref))
            elif state == 3 or state == 4:  # insertion
                self.tindel_lst.append((tpos, ref, alt))
            tpos += ref_len 
            qpos += alt_len


def load_loci_lst(
    region_file: str, 
) -> Tuple[List[str], List[Tuple[str, int, int]]]:

    loci_lst = []
    for line in open(region_file).readlines():
        chrom, start, end = line.strip().split()
        loci_lst.append((chrom, int(start), int(end)))
    return loci_lst


def get_tname2tsize(bam_file: str) -> Tuple[List[str], Dict[str, int]]:

    tname2tsize = {}
    alignments = pysam.AlignmentFile(bam_file, "rb")
    bam_header_lst = str(alignments.header).strip().split("\n")
    for h in bam_header_lst:
        if h.startswith("@SQ"):
            _tag, tname, tsize = h.split("\t")
            tname = tname.replace("SN:", "")
            tsize = tsize.replace("LN:", "")
            tname2tsize[tname] = int(tsize)
    alignments.close()
    return tname2tsize


def dump_ccs_statistics(
    bam_file: str,
    region_file: str,
    min_mapq: int,
    out_file: str
): 

    o = open(out_file, "w")    
    loci_lst = load_loci_lst(region_file)
    alignments = pysam.AlignmentFile(bam_file, "rb")
    for loci in loci_lst:
        for line in alignments.fetch(*loci):
            ccs = BAM(line)
            if not ccs.is_primary:
                continue
            
            if ccs.mapq < min_mapq:
                continue
            o.write(">{}\n{}\n".format(ccs.qname, ccs.qseq)) 
    o.close() 


def main():
    options = parse_args(sys.argv)
    dump_ccs_statistics(
        options.bam, 
        options.region_list,
        options.min_mapq,
        options.out
    )
    sys.exit(0)


if __name__ == "__main__":
    main()

