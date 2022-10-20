#!/usr/bin/env python3

import re
import sys
import pyfastx
import pyabpoa 
import argparse
import multiprocessing as mp
from typing import Dict, List, Tuple, Set

complementary_base_pair = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} 


def parse_args(args):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="mismatch repair candidates",
    )
    parser.add_argument(
        "--ccs",
        type=str,
        required=True,
        help="ccs FASTQ.gz file",
    )
    parser.add_argument(
        "--sscs",
        type=str,
        required=True,
        help="sscs FASTQ.gz file",
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
        "--output",
        type=str,
        required=True,
        help="file to return",
    )
    args = args[1:]
    return parser.parse_args(args)


def load_ccs(infile):

    counter = 0
    qname2dcs = {}
    qname2fwd = {}
    qname2rev = {}
    qname_hsh = {}
    for line in open(infile):
        if line.startswith("qname"):
            continue
        qname, strand, hbit, ref, alt, qbase, qbq, qpos= line.strip().split()
        qbq_lst = qbq.split(",")
        if qbq_lst.count("93") == len(qbq_lst):
            continue
        zmw = qname.split("/")[1]
        dcs = "{}/ccs".format(qname)
        fwd = "{}/ccs/fwd".format(qname)
        rev = "{}/ccs/rev".format(qname)
        qpos_lst = [int(_) for _ in qpos.split(",")]
        qname2dcs[qname] = dcs
        qname2fwd[qname] = fwd
        qname2rev[qname] = rev
        qname_hsh[qname] = [strand, list(hbit), list(ref), list(alt), list(qbase), qbq.split(","), qpos_lst] 

    return qname2dcs, qname2fwd, qname2rev, qname_hsh


def get_compseq(seq: str) -> str:
    compseq = "".join([complementary_base_pair[base] for base in seq])
    return compseq


def get_reverse_compseq(seq: str) -> str:
    reverse_compseq = get_compseq(str(seq)[::-1])
    return reverse_compseq


def get_qpos2mpos(seq: str) -> Dict[int, int]:
    i = 0 
    qpos2mpos = {}
    for j, base in enumerate(seq):
        if base == "-":
            continue
        qpos2mpos[i] = j
        i += 1
    return qpos2mpos

def get_mpos2qpos(seq) -> Dict[int, int]:

    i = 0 
    mpos2qpos = {}
    for j, base in enumerate(seq):
        if base == "-":
            continue
        mpos2qpos[j] = i
        i += 1
    return mpos2qpos

def get_mmr_molecule(
    qname,
    dcs_seq, 
    fwd_seq,
    fwd_qual, 
    rev_seq,
    rev_qual,
    values,
    qname2heteroduplex_lst,
):

    heteroduplex_lst = []
    poa = pyabpoa.msa_aligner()
    strand, hbit_lst, ref_lst, alt_lst, qbase_lst, qbq_lst, qpos_lst = values
    if strand == "+":
        seq_lst = [dcs_seq, fwd_seq, get_reverse_compseq(rev_seq)]
        res=poa.msa(seq_lst, out_cons=False, out_msa=True, out_pog="", incr_fn='') 
    elif strand == "-":
        seq_lst = [get_reverse_compseq(dcs_seq), get_reverse_compseq(fwd_seq), rev_seq]
        res=poa.msa(seq_lst, out_cons=False, out_msa=True, out_pog="", incr_fn='') 
    dcs_msa_seq = res.msa_seq[0]
    fwd_msa_seq = res.msa_seq[1]
    rev_msa_seq = res.msa_seq[2]
    dcs_qpos2mpos = get_qpos2mpos(dcs_msa_seq)
    fwd_mpos2qpos = get_mpos2qpos(fwd_msa_seq)
    rev_mpos2qpos = get_mpos2qpos(rev_msa_seq)
    
    for i, qpos in enumerate(qpos_lst):
        ref = ref_lst[i]
        alt = alt_lst[i]
        qbq = qbq_lst[i]
        qbase = qbase_lst[i]
        mpos = dcs_qpos2mpos[qpos] 
        dca_mbase = dcs_msa_seq[mpos]
        fwd_mbase = fwd_msa_seq[mpos] 
        rev_mbase = rev_msa_seq[mpos]
        if dca_mbase == fwd_mbase == rev_mbase:
            continue
        if (fwd_mbase == ref and rev_mbase == alt) or (fwd_mbase == alt and rev_mbase == ref):
            fwd_bq = ord(fwd_qual[fwd_mpos2qpos[mpos]]) - 33
            rev_bq = ord(rev_qual[rev_mpos2qpos[mpos]]) - 33
            heteroduplex_lst.append([qname, ref, alt, "{}:{}".format(qbase, qbq), "{}:{}".format(fwd_mbase, fwd_bq), "{}:{}".format(rev_mbase, rev_bq)])

    if len(heteroduplex_lst) != 0:
        qname2heteroduplex_lst[qname] = heteroduplex_lst
    
    
def dump_mmr_molecules(
    infile: str,
    ccs_file: str,
    sscs_file: str, 
    threads: int,
    out_file: str,
): 
    ccs_fq = pyfastx.Fastq(ccs_file)
    sscs_fq = pyfastx.Fastq(sscs_file)
    qname2dcs, qname2fwd, qname2rev, qname_hsh = load_ccs(infile)
    
    p = mp.Pool(threads)
    manager = mp.Manager()
    qname2heteroduplex_lst = manager.dict()
    get_mmr_molecule_arg_lst = [
        (
            k,
            ccs_fq[qname2dcs[k]].seq,
            sscs_fq[qname2fwd[k]].seq,
            sscs_fq[qname2fwd[k]].qual,
            sscs_fq[qname2rev[k]].seq,
            sscs_fq[qname2rev[k]].qual,
            v,
            qname2heteroduplex_lst
        )
        for k,v in qname_hsh.items()
    ]
    p.starmap(
        get_mmr_molecule, get_mmr_molecule_arg_lst,
    )
    p.close()
    p.join()

    o = open(out_file, "w")
    o.write("{}\t{}\t{}\t{}\t{}\t{}\n".format("qname", "ref", "alt", "ccs_base", "fwd_base", "rev_base"))
    for _, hduplex_lst in qname2heteroduplex_lst.items():
        for hduplex in hduplex_lst:
            o.write("{}\n".format("\t".join(hduplex))) 
    o.close()



def main():
    options = parse_args(sys.argv)
    dump_mmr_molecules(
        options.input, 
        options.ccs,
        options.sscs,
        options.threads,
        options.output, 
    )
    sys.exit(0)


if __name__ == "__main__":
    main()

