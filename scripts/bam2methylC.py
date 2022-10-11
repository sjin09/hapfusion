
#!/usr/bin/env python3

import re
import os
import sys
import math
import time
import pysam
import tabix
import random
import natsort
import argparse
import statistics
import numpy as np
import multiprocessing as mp
from collections import defaultdict
from typing import Set, Dict, List, Tuple


def parse_args(args):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "-i",
        "--bam",
        type=str,
        required=True,
        help="BAM file to read",
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


class BAM:
    def __init__(self, line):
        self.qname = line.query_name
        self.qseq = line.query_sequence
        self.qlen = len(self.qseq)
        self.bq_int_lst = line.query_qualities


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


def dump_methylated_cytosine(
    bam_file: str,
    out_file: str,
):

    o = open(out_file, "w")
    o.write("{}\t{}\t{}\n".format("qname", "pos", "prob"))
    alignments = pysam.AlignmentFile(bam_file, "rb", check_sq=False)
    for line in alignments:
        ccs = BAM(line)
        if line.has_tag("MM") and line.has_tag("ML"):
            ccs.modified_base_hsh = line.modified_bases
            pos_n_tuple_lst = list(line.modified_bases.values())
            if len(pos_n_tuple_lst) == 0:
                o.write("{}\t{}\t{}\n".format(ccs.qname, ".", "."))
                continue

            p_lst = []
            pos_lst = []
            for (pos, n) in pos_n_tuple_lst[0]:
                p = n/float(256)
                p_lst.append(str(p))
                pos_lst.append(str(pos))
            o.write("{}\t{}\t{}\n".format(ccs.qname, ",".join(pos_lst), ",".join(p_lst)))
    o.close()
            

def main():
    options = parse_args(sys.argv)
    dump_methylated_cytosine(
        options.bam, 
        options.out
    )
    sys.exit(0)


if __name__ == "__main__":
    main()
