import sys
import pysam
import bisect
import pyfastx
import argparse
import himut.bamlib
import himut.caller
import himut.vcflib
from collections import defaultdict
from typing import Set, Dict, List, Tuple
bit_complement_hsh = {"0": "1", "1": "0", "-": "-"}

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
        "--vcf",
        type=str,
        required=True,
        help="phased VCF file to read",
    )    
    parser.add_argument(
        "--bed",
        type=str,
        required=True,
        help="BED file to read",
    )
    parser.add_argument(
        "--sample",
        type=str,
        required=True,
        help="sample id",
    )
    args = args[1:]
    return parser.parse_args(args)


def loci_lst(bed_file):

    chrom_lst = []
    chrom2loci_lst = defaultdict(list)
    for line in open(bed_file):
        chrom, start, end = line.strip().split() 
        loci = (chrom, int(start), int(end)) 
        chrom2loci_lst[chrom].append(loci)
        chrom_lst.append(chrom)
    return chrom_lst, chrom2loci_lst


def get_prdm9_state(
    bam_file: str, 
    loci: Tuple[str, int, int],
    hbit_lst: List[str],
    hpos_lst: List[int],
    hetsnp_lst: List[Tuple[str, int, str, str]]
):

    _chrom, start, end = loci 
    hap2qseq_lst = {"0": [], "1":[], ".": []}
    alignments = pysam.AlignmentFile(bam_file, "rb")
    for i in alignments.fetch(*loci): # iterate through reads
        ccs = himut.bamlib.BAM(i) 
        if start <= ccs.tstart and ccs.tend <= end: # within prdm9
            continue
        elif start <= ccs.tstart and ccs.tend > end: # extends beyond prdm9
            continue
        elif start > ccs.tstart and ccs.tend > end: # spans prdm9 
            ccs_hap = himut.haplib.get_ccs_hap(ccs, hbit_lst, hpos_lst, hetsnp_lst)  
            hap2qseq_lst[ccs_hap].append(ccs)
        elif start > ccs.tstart and ccs.tend <= end: # starts before prdm9 
            continue
    
    if len(hap2qseq_lst["0"]) > 0 and len(hap2qseq_lst["1"]) > 0:
        return True, hap2qseq_lst
    else:
        return False, hap2qseq_lst

def dump_phased_prdm9_sequence(
    sample,
    bam_file,
    vcf_file,
    bed_file,
):

    chrom_lst, chrom2loci_lst = loci_lst(bed_file)
    _, tname2tsize = himut.bamlib.get_tname2tsize(bam_file)
    (
        chrom2ps2hbit_lst,
        chrom2ps2hpos_lst,
        chrom2ps2hetsnp_lst,
        _,
    ) = himut.vcflib.load_phased_hetsnps(vcf_file, chrom_lst, tname2tsize)


    hpos2ps = {}
    hpos_lst = []
    chrom = chrom_lst[0]
    loci = chrom2loci_lst[chrom][0]
    ps2hbit_lst = chrom2ps2hbit_lst[chrom]
    ps2hpos_lst = chrom2ps2hpos_lst[chrom]
    ps2hetsnp_lst = chrom2ps2hetsnp_lst[chrom]
    for ps, ps_hpos_lst in ps2hpos_lst.items():
       for ps_hpos in ps_hpos_lst:
            hpos2ps[ps_hpos] = ps  
            hpos_lst.append(ps_hpos)
 
    hap_lst = ["0", "1"]
    chrom, start, end = loci 
    hap2ccs_lst = defaultdict(list)
    idx = bisect.bisect_right(hpos_lst, start)
    jdx = bisect.bisect_right(hpos_lst, end)
    hpos_subset_lst = [hpos_lst[kdx] for kdx in range(idx, jdx)]
    phase_set_lst = list(set([hpos2ps[hpos] for hpos in hpos_subset_lst]))
    if len(phase_set_lst) == 1:
        ps = phase_set_lst[0]        
        hbit_lst = ps2hbit_lst[ps]
        hpos_lst = ps2hpos_lst[ps]
        hetsnp_lst = ps2hetsnp_lst[ps]
        alignments = pysam.AlignmentFile(bam_file, "rb")
        for i in alignments.fetch(*loci): # iterate through reads
            ccs = himut.bamlib.BAM(i) 
            ccs.hbit = ""
            ccs.cs2tpos2qbase()
            idx = bisect.bisect_right(hpos_lst, ccs.tstart)
            jdx = bisect.bisect_right(hpos_lst, ccs.tend)
            hetsnp_subset_lst = [hetsnp_lst[kdx] for kdx in range(idx, jdx)]
            h0_hbit = "".join([hbit_lst[kdx] for kdx in range(idx, jdx)])
            h1_hbit = "".join([bit_complement_hsh[hbit] for hbit in h0_hbit])
            for hetsnp in hetsnp_subset_lst:  # get read haplotype bits
                qbase = ccs.tpos2qbase[hetsnp[0]][0]
                if qbase == hetsnp[1]:  # ref
                    ccs.hbit += "0"
                elif qbase == hetsnp[2]:  # alt
                    ccs.hbit += "1"
                else:
                    ccs.hbit += "-"

            if h0_hbit == ccs.hbit:
                ccs.hap = "0"
            elif h1_hbit == ccs.hbit:
                ccs.hap = "1"
            else:
                ccs.hap = "."
            hap2ccs_lst[ccs.hap].append(ccs)
        
        for hap in hap_lst:
            tstart_lst = []
            tend_lst = []
            ccs_lst = hap2ccs_lst[hap] 
            for ccs in ccs_lst:
                tend_lst.append(ccs.tend)
                tstart_lst.append(ccs.tstart)
            
            tstart = min(tstart_lst)
            tend = max(tend_lst)
            if tstart < start and end < tend:
                o = open("{}.h{}.prdm9.ccs.fastq".format(sample, hap), "w")
                for ccs in ccs_lst:
                    qbq = "".join([chr(bq+33) for bq in ccs.bq_int_lst])
                    o.write("@{}\n{}\n+\n{}\n".format(ccs.qname, ccs.qseq, qbq)) 
                o.close()
            print(hap, tstart, tend, tend-tstart)

def main():
    options = parse_args(sys.argv)
    dump_phased_prdm9_sequence(
        options.sample, 
        options.bam,
        options.vcf,
        options.bed,
    )
    sys.exit(0)


if __name__ == "__main__":
    main()

