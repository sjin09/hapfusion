
#!/usr/bin/env python3

import re
import sys
import pysam
import tabix
import cyvcf2
import bisect
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
        help="BAM file",
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
        "--vcf",
        type=str,
        required=True,
        help="phased deepvariant VCF file",
    )
    parser.add_argument(
        "--min_mapq",
        type=int,
        default=60,
        required=False,
        help="minimum mapping quality score",
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
        ## self.qlen = len(self.qseq)
        self.mapq = line.mapping_quality
        self.bq_int_lst = line.query_qualities
        ## self.qv = sum(self.bq_int_lst)/self.qlen
        self.strand = "+" if line.is_forward else "-"
        ## self.hbq_proportion = self.bq_int_lst.count(93)/float(self.qlen)
        self.query_alignment_length = self.qend - self.qstart
        ## self.query_alignment_proportion = self.query_alignment_length/float(self.qlen)
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


    def cs2tpos2qbase(self):
        self.cs2tuple()
        tpos = self.tstart
        qpos = self.qstart
        self.tpos2qbase = {}
        for cstuple in self.cstuple_lst:
            state, ref, alt, ref_len, alt_len, = cstuple
            if state == 1:  # match
                for i, alt_base in enumerate(alt):
                    self.tpos2qbase[tpos + i + 1] = (qpos+i, alt_base, self.bq_int_lst[qpos + i])
            elif state == 2:  # sub
                self.tpos2qbase[tpos + 1] = (qpos, alt, self.bq_int_lst[qpos])
            elif state == 3:  # insertion
                pass
            elif state == 4:  # deletion
                for j in range(len(ref[1:])):
                    self.tpos2qbase[tpos + j + 1] = (qpos, "-", 0)
            tpos += ref_len
            qpos += alt_len
        return self.tpos2qbase



class VCF:
    def __init__(self, line):
        arr = line.strip().split()
        self.chrom = arr[0]
        self.pos = int(arr[1])
        self.id = arr[2]
        self.ref = arr[3]
        self.alt_lst = arr[4].split(",")
        self.qual = arr[5]
        self.qual = float(self.qual) if self.qual != "." else self.qual
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


def load_haplotype_block(
    vcf_file: str, 
    chrom: str,
    chrom_len: int,
) -> Dict[str, List[List[Tuple[int, str]]]]:

    hidx = 0  
    hidx2hetsnp = {}
    hetsnp2hidx = {}
    hblock_lst = defaultdict(list)
    phase_set2hblock = defaultdict(list) 
    if vcf_file.endswith(".vcf"):
        for line in open(vcf_file):
            if line.startswith("#"):
                continue
            v = VCF(line)
            if chrom == v.chrom: 
                if v.sample_gt == "0|1" or v.sample_gt == "1|0":
                    hetsnp = (v.pos, v.ref, v.alt)
                    hidx2hetsnp[hidx] = hetsnp
                    hetsnp2hidx[hetsnp] = hidx
                    hstate = v.sample_gt.split("|")[0]                    
                    phase_set2hblock[v.sample_phase_set].append((hidx, hstate))
                hidx += 1 
    elif vcf_file.endswith(".bgz"):
        tb = tabix.open(vcf_file)
        records = tb.query(chrom, 0, chrom_len)
        for hidx, record in enumerate(records):
            v = VCF("\t".join(record))           
            if v.sample_gt == "0|1" or v.sample_gt == "1|0":
                hetsnp = (v.pos, v.ref, v.alt)
                hidx2hetsnp[hidx] = hetsnp
                hetsnp2hidx[hetsnp] = hidx                
                hstate = v.sample_gt.split("|")[0]
                phase_set2hblock[v.sample_phase_set].append((hidx, hstate))
    hblock_lst = [hblock for hblock in phase_set2hblock.values()]
    return hblock_lst, hidx2hetsnp, hetsnp2hidx


def get_haplotype_block(
    chrom: str,
    chrom_len: int, 
    vcf_file: str, 
):
    filtered_hblock_lst = []
    hblock_lst, hidx2hetsnp, _ = load_haplotype_block(vcf_file, chrom, chrom_len)
    for hblock in hblock_lst:
        ipos = 0
        filtered_hblock = []
        for (hidx, hstate) in hblock:
            jpos = hidx2hetsnp[hidx][0]
            if jpos - ipos > 1:
                filtered_hblock.append((hidx, hstate))
            ipos = jpos
        if len(filtered_hblock) > 1:
            filtered_hblock_lst.append(filtered_hblock)

    bidx2loci = {}
    bidx2hbit = {}
    bidx2hetsnp_lst = defaultdict(list)
    for bidx, hblock in enumerate(filtered_hblock_lst):
        hstate_lst = []
        for (hidx, hstate) in hblock:
            hstate_lst.append(hstate) 
            hetsnp = hidx2hetsnp[hidx]
            bidx2hetsnp_lst[bidx].append(hetsnp)
        hbit = "".join(hstate_lst) 
        bstart = bidx2hetsnp_lst[bidx][0][0]
        bend = bidx2hetsnp_lst[bidx][-1][0]
        blength = bend - bstart 
        if blength > 1000:
            bidx2hbit[bidx] = hbit 
            bidx2loci[bidx] = (chrom, bstart, bend)
    return bidx2loci, bidx2hbit, bidx2hetsnp_lst


def get_mmr_candidates(
    chrom: str,
    chrom_len: int,
    bam_file: str,
    vcf_file: str,
    min_mapq: int,
    chrom2ccs_lst: Dict[str, str] 
):
  
    ccs_lst = [] 
    alignments = pysam.AlignmentFile(bam_file, "rb")
    bidx2loci, bidx2hbit, bidx2hetsnp_lst = get_haplotype_block(chrom, chrom_len, vcf_file)
    for bidx, loci in bidx2loci.items():
        block_hbit = bidx2hbit[bidx]
        hetsnp_lst = bidx2hetsnp_lst[bidx]
        hpos_lst = [hetsnp[0] for hetsnp in hetsnp_lst]
        for line in alignments.fetch(*loci):
            ccs = BAM(line)
            if not ccs.is_primary:
                continue
            
            if ccs.mapq < min_mapq:
                continue
            
            idx = bisect.bisect_right(hpos_lst, ccs.tstart)
            jdx = bisect.bisect_right(hpos_lst, ccs.tend)
            # block_subset_hbit = block_hbit[idx:jdx]
            hetsnp_subset_lst = hetsnp_lst[idx:jdx]
            hetsnp_subset_cnt = len(hetsnp_subset_lst)
            if hetsnp_subset_cnt < 3:
                continue
            
            qbq_lst = [] 
            ref_lst = []
            alt_lst = []
            qpos_lst = []
            qbase_lst = []
            ccs_hbit = ""
            ccs.cs2tpos2qbase()
            for (tpos, ref, alt)  in hetsnp_subset_lst:
                qpos, qbase, qbq = ccs.tpos2qbase[tpos]
                if qbase == ref:
                    ccs_hbit += "0"
                    ref_lst.append(ref)
                    alt_lst.append(alt)
                    qbq_lst.append(str(qbq))
                    qpos_lst.append(str(qpos))
                    qbase_lst.append(qbase)
                elif qbase == alt:
                    ccs_hbit += "1"
                    ref_lst.append(ref)
                    alt_lst.append(alt)
                    qbq_lst.append(str(qbq))
                    qpos_lst.append(str(qpos))
                    qbase_lst.append(qbase)
            ccs_lst.append(
                [
                    ccs.qname, 
                    ccs.strand,
                    ccs_hbit,
                    "".join(ref_lst), 
                    "".join(alt_lst), 
                    "".join(qbase_lst),
                    ",".join(qbq_lst),
                    ",".join(qpos_lst), 
                ]
            ) 
    chrom2ccs_lst[chrom] = ccs_lst
    

def dump_mmr_candidates(
    bam_file: str,
    vcf_file: str, 
    region: str,
    region_file: str,
    min_mapq: int,
    threads: int,
    out_file: str
): 
    
    p = mp.Pool(threads)
    manager = mp.Manager()
    chrom2ccs_lst = manager.dict()
    tname2tsize = get_tname2tsize(bam_file)
    chrom_lst = load_chrom_lst(region, region_file)
    get_mmr_candidates_arg_lst = [
        (
            chrom,
            tname2tsize[chrom],
            bam_file,
            vcf_file,
            min_mapq,
            chrom2ccs_lst
        )
        for chrom in chrom_lst
    ]
    p.starmap(
        get_mmr_candidates, get_mmr_candidates_arg_lst,
    )
    p.close()
    p.join()
  
    o = open(out_file, "w")
    o.write("{}\n".format("\t".join(["qname", "strand", "hbit", "ref", "alt", "qbase", "qbq", "qpos"])))
    for chrom in chrom_lst:
        for (qname, strand, hbit, ref, alt, qbase, qbq, qpos) in chrom2ccs_lst[chrom]: 
            o.write("{}\n".format("\t".join([qname, strand, hbit, ref, alt, qbase, qbq, qpos])))
    o.close() 


def main():
    options = parse_args(sys.argv)
    dump_mmr_candidates(
        options.bam, 
        options.vcf,
        options.region, 
        options.region_list,
        options.min_mapq,
        options.threads, 
        options.out
    )
    sys.exit(0)


if __name__ == "__main__":
    main()
