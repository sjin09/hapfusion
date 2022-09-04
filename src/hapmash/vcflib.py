import os
import json
import tabix
import cyvcf2
import natsort
import numpy as np 
import hapmash.util
import hapmash.bamlib
import hapmash.phaselib
from datetime import datetime
from collections import defaultdict
from typing import Dict, List, Set, Tuple


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


def get_sample(vcf_file: str) -> str:
    if vcf_file.endswith(".vcf"):
        for line in open(vcf_file).readlines():
            if line.startswith("#CHROM"):
                sample = line.strip().split()[-1]
                return sample
    elif vcf_file.endswith(".bgz"):
        v = cyvcf2.VCF(vcf_file)
        sample = v.samples[0]
        return sample


def get_hapmash_vcf_header(
    version: str,
    bam_file: str,
) -> str:

    vcf_header_lst = [
        "##fileformat=VCFv4.2",
        '##FILTER=<ID=PASS,Description="All filters passed">',
        "##fileDate={}".format(datetime.now().strftime("%d%m%Y")),
        "##source=hapmash",
        "##source_version={}".format(version),
        '##content=hapmash gene conversions',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '##FORMAT=<ID=BQ,Number=1,Type=Integer,Description="Average base quality">',
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">',
        '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Read depth for each allele">',
        '##FORMAT=<ID=VAF,Number=A,Type=Float,Description="Variant allele fractions">',
        '##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">',
    ]
    tname2tsize = hapmash.bamlib.get_tname2tsize(bam_file)
    tname_lst = natsort.natsorted(tname2tsize.keys())
    for tname in tname_lst:
        vcf_header_lst.append("##contig=<ID={},length={}>".format(tname, tname2tsize[tname]))
    vcf_header_lst.append(
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}".format(hapmash.bamlib.get_sample(bam_file))
    )
    vcf_header = "\n".join(vcf_header_lst)
    return vcf_header


def get_phased_vcf_header(
    bam_file: str,
    version: str,
) -> str:

    vcf_header_lst = [
        "##fileformat=VCFv4.2",
        '##FILTER=<ID=PASS,Description="All filters passed">',
        "##fileDate={}".format(datetime.now().strftime("%d%m%Y")),
        "##source=hapmash",
        "##source_version={}".format(version),
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Conditional genotype quality">',
        '##FORMAT=<ID=BQ,Number=1,Type=Integer,Description="Average base quality">',
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">',
        '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Read depth for each allele">',
        '##FORMAT=<ID=VAF,Number=A,Type=Float,Description="Variant allele fractions">',
        '##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred-scaled genotype likelihoods rounded to the closest integer">',
        '##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">',
    ]
    tname_lst, tname2tsize = hapmash.bamlib.get_tname2tsize(bam_file)
    for tname in tname_lst:
        vcf_header_lst.append("##contig=<ID={},length={}>".format(tname, tname2tsize[tname]))

    vcf_header_lst.append(
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}".format(hapmash.bamlib.get_sample(bam_file))
    )
    vcf_header = "\n".join(vcf_header_lst)
    return vcf_header


def load_snp_indel(
    chrom: str,
    vcf_file: str
) -> Set[Tuple[str, int, str, str]]:

    homsnp_set = set()
    homindel_set = set()
    hetindel_set = set()
    for line in open(vcf_file):
        if line.startswith("#"):
            continue
        v = VCF(line)
        if v.chrom != chrom:
            continue
        if v.is_pass and v.is_snp and v.is_biallelic and v.sample_gt == "1/1":
            homsnp_set.add((v.chrom, v.pos, v.ref, v.alt))
        elif v.is_pass and v.is_indel and v.is_biallelic and v.sample_gt == "0/1":
            hetindel_set.add((v.chrom, v.pos, v.ref, v.alt))
        elif v.is_pass and v.is_indel and v.is_biallelic and v.sample_gt == "1/1":
            homindel_set.add((v.chrom, v.pos, v.ref, v.alt))
    return homsnp_set, hetindel_set, homindel_set
      
                     
def load_bgz_snp_indel(
    loci: Tuple[str, int, int],
    vcf_file: str
) -> Set[Tuple[str, int, str, str]]:

    het_set = set()
    hom_set = set()
    tb = tabix.open(vcf_file)
    records = tb.query(*loci)
    for record in records:
        v = VCF("\t".join(record))            
        # if v.is_pass and v.is_snp and v.is_biallelic and v.sample_gt == "0/1":
        #     het_set.add((v.chrom, v.pos, v.ref, v.alt))
        # elif v.is_pass and v.is_snp and v.is_biallelic and v.sample_gt == "1/1":
        #     hom_set.add((v.chrom, v.pos, v.ref, v.alt))
        # elif v.is_pass and v.is_indel and v.is_biallelic and v.sample_gt == "0/1":
        #     het_set.add((v.chrom, v.pos, v.ref, v.alt))
        # elif v.is_pass and v.is_indel and v.is_biallelic and v.sample_gt == "1/1":
        #     hom_set.add((v.chrom, v.pos, v.ref, v.alt))
        if v.is_pass and v.is_biallelic and v.sample_gt == "0/1":
            het_set.add((v.chrom, v.pos, v.ref, v.alt))
        elif v.is_pass and v.is_biallelic and v.sample_gt == "1/1":
            hom_set.add((v.chrom, v.pos, v.ref, v.alt))
    return het_set, hom_set, 


def load_sbs(
    chrom: str,
    vcf_file: str
) -> Set[Tuple[str, int, str, str]]:

    sbs_set = set()
    for line in open(vcf_file):
        if line.startswith("#"):
            continue
        v = VCF(line)
        if v.chrom != chrom:
            continue
        if v.is_pass and v.is_snp and v.is_biallelic:
            sbs_set.add((v.chrom, v.pos, v.ref, v.alt))
    return sbs_set 
      
                     
def load_bgz_sbs(
    loci: Tuple[str, int, int],
    vcf_file: str
) -> Set[Tuple[str, int, str, str]]:

    sbs_set = set()
    tb = tabix.open(vcf_file)
    records = tb.query(*loci)
    for record in records:
        v = VCF("\t".join(record))            
        if v.is_pass and v.is_snp and v.is_biallelic:
            sbs_set.add((v.chrom, v.pos, v.ref, v.alt))
    return sbs_set


def load_hetsnps(
    vcf_file: str,
    chrom: str,
    chrom_len: int
) -> Dict[str, List[Tuple[str, int, str, str]]]:

    hetsnp_lst = []
    hidx2hetsnp = {}
    hetsnp2hidx = {}
    if vcf_file.endswith(".vcf"):
        for line in open(vcf_file):
            if line.startswith("#"):
                continue
            v = VCF(line)
            if v.is_snp and v.is_pass and v.is_biallelic and (v.sample_gt == "0/1" or v.sample_gt == "1/0"):
                hetsnp_lst.append((v.chrom, v.pos, v.ref, v.alt))
    elif vcf_file.endswith(".bgz"):
        tb = tabix.open(vcf_file)
        records = tb.query(chrom, 0, chrom_len)
        for record in records:
            v = VCF("\t".join(record))            
            if v.is_snp and v.is_pass and v.is_biallelic and (v.sample_gt == "0/1" or v.sample_gt == "1/0"):
                hetsnp_lst.append((v.chrom, v.pos, v.ref, v.alt))
    for hidx, hetsnp in enumerate(hetsnp_lst):
        hidx2hetsnp[hidx] = hetsnp
        hetsnp2hidx[hetsnp] = hidx
    return hetsnp_lst, hidx2hetsnp, hetsnp2hidx


def load_haplotype_block(
    vcf_file: str, 
    chrom: str,
    chrom_len: int,
) -> Dict[str, List[List[Tuple[int, str]]]]:

    hidx = 0  
    state = 0
    hidx2hetsnp = {}
    hetsnp2hidx = {}
    hblock_lst = defaultdict(list)
    phase_set2hblock = defaultdict(list) 
    if vcf_file.endswith(".vcf"):
        for line in open(vcf_file):
            if line.startswith("#"):
                continue
            v = VCF(line)
            if chrom != v.chrom:
                if state == 0:
                    continue
                elif state == 1:
                    break
            elif chrom == v.chrom: 
                if state == 0:
                    state = 1
                if v.is_pass:
                    if v.is_snp and v.is_biallelic:
                        snp = (v.chrom, v.pos, v.ref, v.alt)
                        if v.sample_gt == "0|1" or v.sample_gt == "1|0":
                            hidx2hetsnp[hidx] = snp
                            hetsnp2hidx[snp] = hidx
                            hstate = v.sample_gt.split("|")[0]                    
                            phase_set2hblock[v.sample_phase_set].append((hidx, hstate))
                            hidx += 1 
    elif vcf_file.endswith(".bgz"):
        tb = tabix.open(vcf_file)
        records = tb.query(chrom, 0, chrom_len)
        for record in records:
            v = VCF("\t".join(record)) 
            if v.is_pass:
                if v.is_snp and v.is_biallelic:
                    snp = (v.chrom, v.pos, v.ref, v.alt)
                    if v.sample_gt == "0|1" or v.sample_gt == "1|0":
                        hidx2hetsnp[hidx] = snp
                        hetsnp2hidx[snp] = hidx
                        hstate = v.sample_gt.split("|")[0]                    
                        phase_set2hblock[v.sample_phase_set].append((hidx, hstate))
                        hidx += 1 
    hblock_lst = [hblock for hblock in phase_set2hblock.values()]
    return hblock_lst, hidx2hetsnp, hetsnp2hidx


def get_phased_hetsnps(
    vcf_file: str, 
    chrom: str,
    chrom_len: int, 
):
    
    hidx = 0  
    state = 0
    hidx2hetsnp = {}
    hetsnp2hidx = {}
    hblock_lst = defaultdict(list)
    phase_set2hblock = defaultdict(list) 
    if vcf_file.endswith(".vcf"):
        for line in open(vcf_file):
            if line.startswith("#"):
                continue
            v = VCF(line)
            if chrom != v.chrom:
                if state == 0:
                    continue
                elif state == 1:
                    break
            elif chrom == v.chrom: 
                if state == 0:
                    state = 1
                if v.is_pass:
                    if v.is_snp and v.is_biallelic:
                        snp = (v.chrom, v.pos, v.ref, v.alt)
                        if v.sample_gt == "0|1" or v.sample_gt == "1|0":
                            hidx2hetsnp[hidx] = snp
                            hetsnp2hidx[snp] = hidx
                            hstate = v.sample_gt.split("|")[0]                    
                            phase_set2hblock[v.sample_phase_set].append((hidx, hstate))
                            hidx += 1 
    elif vcf_file.endswith(".bgz"):
        tb = tabix.open(vcf_file)
        records = tb.query(chrom, 0, chrom_len)
        for record in records:
            v = VCF("\t".join(record)) 
            if v.is_pass:
                if v.is_snp and v.is_biallelic:
                    snp = (v.chrom, v.pos, v.ref, v.alt)
                    if v.sample_gt == "0|1" or v.sample_gt == "1|0":
                        hidx2hetsnp[hidx] = snp
                        hetsnp2hidx[snp] = hidx
                        hstate = v.sample_gt.split("|")[0]                    
                        phase_set2hblock[v.sample_phase_set].append((hidx, hstate))
                        hidx += 1 
    hblock_lst = [hblock for hblock in phase_set2hblock.values()]

    hetsnp2bidx = {}
    hidx2hstate = {} 
    filtered_hblock_lst = []
    for hblock in hblock_lst:
        ipos = 0
        filtered_hblock = []
        for (hidx, hstate) in hblock:
            jpos = hidx2hetsnp[hidx][1]
            if jpos - ipos > 1:
                filtered_hblock.append((hidx, hstate))
            ipos = jpos
        if len(filtered_hblock) > 1:
            filtered_hblock_lst.append(filtered_hblock)
   
    hidx2bidx = {}
    hetsnp2hstate = {}
    phased_hetsnp_lst = []
    for bidx, hblock in enumerate(filtered_hblock_lst):
        for (hidx, hstate) in hblock:
            hidx2bidx[hidx] = bidx
            hidx2hstate[hidx] = hstate 
            hetsnp = hidx2hetsnp[hidx]
            hetsnp2bidx[hetsnp] = bidx 
            hetsnp2hstate[hetsnp] = hstate
            phased_hetsnp_lst.append(hetsnp)
    phased_hetsnp_set = set(natsort.natsorted(phased_hetsnp_lst))
    phased_hpos_lst = [hetsnp[1] for hetsnp in phased_hetsnp_lst]
    return phased_hpos_lst, phased_hetsnp_lst, phased_hetsnp_set, filtered_hblock_lst, hidx2bidx, hidx2hstate, hidx2hetsnp, hetsnp2bidx, hetsnp2hidx, hetsnp2hstate


def dump_phased_hetsnps(
    bam_file: str,
    vcf_file: str,
    min_qual: int,
    chrom_lst: List[str],
    tname2tsize: Dict[str, int],
    chrom2hblock_lst: List[List[Tuple[int, int]]],
    version: str,
    out_file: str,
) -> None:

    hetsnp2hstate = {}
    hetsnp2phase_set = {}
    chrom2hetsnp_lst = defaultdict(list)
    chrom2hidx2hetsnp = defaultdict(dict)
    chrom2hetsnp2hidx = defaultdict(dict)
    if vcf_file.endswith(".vcf"): # load germline mutations
        for line in open(vcf_file):
            if line.startswith("#"):
                continue
            v = VCF(line)
            if v.is_snp and v.is_pass and v.is_biallelic and (v.sample_gt == "0/1" or v.sample_gt == "1/0"):
                chrom2hetsnp_lst[v.chrom].append((v.chrom, v.pos, v.ref, v.alt))
    elif vcf_file.endswith(".bgz"):
        tb = tabix.open(vcf_file)
        for chrom in chrom_lst:
            records = tb.query(chrom, 0, tname2tsize[chrom])
            for record in records:
                v = VCF("\t".join(record))            
                if v.is_snp and v.is_pass and v.is_biallelic and (v.sample_gt == "0/1" or v.sample_gt == "1/0"):
                    chrom2hetsnp_lst[v.chrom].append((v.chrom, v.pos, v.ref, v.alt))
    for chrom, hetsnp_lst in chrom2hetsnp_lst.items():
        for hidx, hetsnp in enumerate(hetsnp_lst):
            chrom2hidx2hetsnp[chrom][hidx] = hetsnp
            chrom2hetsnp2hidx[chrom][hetsnp] = hidx
   
    for chrom, hblock_lst in chrom2hblock_lst.items():
        for hblock in hblock_lst:
            phase_set = chrom2hidx2hetsnp[chrom][hblock[0][0]][1]
            for (hidx, hstate) in hblock:
                phased_hetsnp = chrom2hidx2hetsnp[chrom][hidx]
                hetsnp2hstate[phased_hetsnp] = hstate
                hetsnp2phase_set[phased_hetsnp] = phase_set 
    
    o = open(out_file, "w") # return phased hetsnps and germline mutations
    o.write("{}\n".format(hapmash.vcflib.get_phased_vcf_header(bam_file, version))) 
    if vcf_file.endswith(".vcf"):
        for line in open(vcf_file):
            if line.startswith("#"):
                continue
            v = VCF(line)
            if v.is_pass:
                fmt, sample_fmt = line.strip().split()[-2:]
                if v.is_snp and v.is_biallelic and (v.sample_gt == "0/1" or v.sample_gt == "1/0"):
                    hetsnp = (v.chrom, v.pos, v.ref, v.alt)
                    if v.qual > min_qual:
                        if hetsnp in hetsnp2hstate:
                            phase_set = hetsnp2phase_set[hetsnp]
                            sample_gt = "0|1" if hetsnp2hstate[hetsnp] == "0" else "1|0"
                            sample_fmt_arr = sample_fmt.split(":")
                            sample_fmt_arr[0] = sample_gt
                            o.write("{}\t{}\t.\t{}\t{}\t{}\tPASS\t.\t{}:PS\t{}:{}\n".format(v.chrom, v.pos, v.ref, v.alt, v.qual, fmt, ":".join(sample_fmt_arr), phase_set))
                        else:
                            o.write("{}\t{}\t.\t{}\t{}\t{}\tPASS\t.\t{}:PS\t{}:.\n".format(v.chrom, v.pos, v.ref, v.alt, v.qual, fmt, sample_fmt))
                    else:
                        o.write("{}\t{}\t.\t{}\t{}\t{}\tPASS\t.\t{}:PS\t{}:.\n".format(v.chrom, v.pos, v.ref, v.alt, v.qual, fmt, sample_fmt))
                else:
                    if v.is_biallelic:
                        o.write("{}\t{}\t.\t{}\t{}\t{}\tPASS\t.\t{}:PS\t{}:.\n".format(v.chrom, v.pos, v.ref, v.alt, v.qual, fmt, sample_fmt))
                    else: 
                        o.write("{}\t{}\t.\t{}\t{}\t{}\tPASS\t.\t{}:PS\t{}:.\n".format(v.chrom, v.pos, v.ref, ",".join(v.alt_lst), v.qual, fmt, sample_fmt))
    elif vcf_file.endswith(".bgz"):
        tb = tabix.open(vcf_file)
        for chrom in chrom_lst:
            records = tb.query(chrom, 0, tname2tsize[chrom])
            for record in records:
                line = "\t".join(record)
                v = VCF(line)            
                fmt, sample_fmt = line.strip().split()[-2:]
                if v.is_pass:
                    if v.is_snp and v.is_biallelic and (v.sample_gt == "0/1" or v.sample_gt == "1/0"):
                        hetsnp = (v.chrom, v.pos, v.ref, v.alt)
                        if v.qual > min_qual:
                            if hetsnp in hetsnp2hstate:
                                phase_set = hetsnp2phase_set[hetsnp]
                                sample_gt = "0|1" if hetsnp2hstate[hetsnp] == "0" else "1|0"
                                sample_fmt_arr = sample_fmt.split(":")
                                sample_fmt_arr[0] = sample_gt
                                o.write("{}\t{}\t.\t{}\t{}\t{}\tPASS\t.\t{}:PS\t{}:{}\n".format(v.chrom, v.pos, v.ref, v.alt, v.qual, fmt, ":".join(sample_fmt_arr), phase_set))
                            else:
                                o.write("{}\t{}\t.\t{}\t{}\t{}\tPASS\t.\t{}:PS\t{}:.\n".format(v.chrom, v.pos, v.ref, v.alt, v.qual, fmt, sample_fmt))
                        else:
                            o.write("{}\t{}\t.\t{}\t{}\t{}\tPASS\t.\t{}:PS\t{}:.\n".format(v.chrom, v.pos, v.ref, v.alt, v.qual, fmt, sample_fmt))
                    else:
                        if v.is_biallelic:
                            o.write("{}\t{}\t.\t{}\t{}\t{}\tPASS\t.\t{}:PS\t{}:.\n".format(v.chrom, v.pos, v.ref, v.alt, v.qual, fmt, sample_fmt))
                        else: 
                            o.write("{}\t{}\t.\t{}\t{}\t{}\tPASS\t.\t{}:PS\t{}:.\n".format(v.chrom, v.pos, v.ref, ",".join(v.alt_lst), v.qual, fmt, sample_fmt))
    o.close()


def dump_hblock_statistics(
    chrom_lst: List[str],
    chrom2hblock_statistics: Dict[str, Tuple[int, int, int, int, int]],
    log_file: str, 
) -> None:

    o = open(log_file, "w")
    o.write(
        "chrom\thetsnp_count\tblock_count\tsmallest_block\tlargest_block\tshortest_block (bp)\tlongest_block (bp)\n"
    )
    for chrom in chrom_lst:
        statistics = chrom2hblock_statistics[chrom]
        o.write(
            "{0}\t{1[0]}\t{1[1]}\t{1[2]}\t{1[3]:,}\t{1[4]:,}\t{1[5]}\n".format(
                chrom, statistics
            )
        )
    o.close()


def dump_recombinantion(
    chrom_lst: List[str], 
    chrom2hapmash_lst: List[str],
    out_file: str
):
    o = open(out_file, "w")
    for chrom in chrom_lst:
        for (ccs, ccs_haplotype, hapmash_hetsnp_lst, wmark_hetsnp_lst, wmark_lst, h0_hbit, h1_hbit, read_hbit) in chrom2hapmash_lst[chrom]:
            o.write("{}:{}\t{}\t{}\n\n{}\n{}\n{}\n{}\n\n".format(ccs, ccs_haplotype, ",".join(hapmash_hetsnp_lst), wmark_hetsnp_lst, wmark_lst, h0_hbit, h1_hbit, read_hbit))
    o.close()


def dump_hetsnp_statistics(
    chrom_lst: List[str], 
    chrom2hetsnp_statistics: List[str],
    out_file: str
):
    o = open(out_file, "w")
    hsh = defaultdict(lambda: 0)
    for chrom in chrom_lst:
        for hetsnp_cnt, read_cnt in chrom2hetsnp_statistics[chrom].items():
            hsh[hetsnp_cnt] += read_cnt
    for hetsnp_cnt in natsort.natsorted(list(hsh.keys())):
            o.write("{}\t{}\n".format(hetsnp_cnt, hsh[hetsnp_cnt]))
    o.close()
                

def dump_recombinantion_statistics(
    chrom_lst: List[str], 
    chrom2hapmash_statistics: Dict[str, List[int]],
    out_file: str
):
    
    row_names = [
    "ccs",
    "lq_ccs",
    "hq_ccs",
    "unphaseable_hq_ccs",
    "phaseable_hq_ccs",
    "unphased_hq_ccs",
    "phased_hq_ccs",
    "hap_consistent_hq_ccs", 
    "hap_inconsistent_hq_ccs",
    "recombinant_ccs"
    ]
    ncol = len(chrom_lst)
    nrow = len(row_names)
    dt = np.zeros((nrow, ncol))
    for idx, chrom in enumerate(chrom_lst):
        for jdx, count in enumerate(chrom2hapmash_statistics[chrom]): 
            dt[jdx][idx] = count
    
    o = open(out_file, "w")
    genome_lst = chrom_lst + ["total"] 
    o.write("{:30}{}\n".format("", "\t".join(genome_lst)))
    for kdx in range(nrow):
        rsum =  str(int(np.sum(dt[kdx])))
        stats = "\t".join([str(int(stat)) for stat in dt[kdx].tolist()] + [rsum])
        o.write("{:30}{}\n".format(row_names[kdx], stats))
    o.close()


      
