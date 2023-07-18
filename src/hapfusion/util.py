import os
import sys
import pysam
import tabix
import cyvcf2
from typing import Dict, List
from collections import defaultdict
base_lst = list("ATGC")
base_set = set(base_lst)
allele_lst = list("ATGC+-")
base2idx = {base: idx for idx, base in enumerate(base_lst)}
idx2base = {idx: base for idx, base in enumerate(base_lst)}
allele2idx = {allele: idx for idx, allele in enumerate(allele_lst)}
idx2allele = {idx: allele for idx, allele in enumerate(allele_lst)}


class FUSION:
    def __init__(self, line):
        arr = line.strip().split()
        coord = arr[0] 
        self.tname = coord.split(":")[0]
        self.tstart, self.tend = coord.split(":")[1].split("-")
        self.tstart = int(self.tstart)
        self.tend = int(self.tend)
        self.qname = arr[1]
        self.phase_set = arr[2]
        self.is_pass = True if arr[3] == "PASS" else False
        self.event = arr[4] 
        self.haplotype = arr[5]
        self.hd = float(arr[6]) if arr[6] != "." else "."
        self.recomb_length = int(arr[7]) if arr[7] != "." else "."
        self.ccs_hbit = arr[8]
        self.h0_hbit = arr[9]
        self.h1_hbit = arr[10]
        self.hetsnps = arr[11]
        # self.denovo_mutations = arr[12]


def exit():
    print("exiting hapfusion")
    sys.exit(0)


def is_bam_file_corrupt(bam_file: str, chrom_lst: List[str]):
    hsh = {}
    alignments = pysam.AlignmentFile(bam_file, "rb")
    for chrom in chrom_lst:
        hsh[chrom] = alignments.count(chrom)

    state = 0
    for chrom in chrom_lst:
        if hsh[chrom] == 0:
            print(
                "bam_file: {} does not have read alignments for chrom: {}".format(
                    bam_file, chrom
                )
            )
            state = 1

    if state == 1:
        print("{} might be corrupted".format(bam_file))
        return 1
    else:
        return 0


def check_bam_file(
    bam_file: str,
    chrom_lst: List[str],
):

    if bam_file is None:
        print("Please provide the path to the BAM file")
        return 1
    else:
        if bam_file.endswith(".bam"):
            idxfile = "{}.bai".format(bam_file)
            if os.path.exists(bam_file) and os.path.exists(idxfile):
                if os.path.getsize(bam_file) != 0 and os.path.getsize(idxfile) != 0:
                    if is_bam_file_corrupt(bam_file, chrom_lst):
                        return 0
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
            return 1


def is_vcf_file_corrupt(
    vcf_file: str, chrom_lst: List[str], tname2tsize: Dict[str, int]
) -> bool:

    hsh = defaultdict(lambda: 0)
    if vcf_file.endswith(".vcf"):
        for line in open(vcf_file):
            if line.startswith("#"):
                continue
            arr = line.strip().split()
            chrom = arr[0]
            ref = arr[3]
            alt_lst = arr[4].split(",")
            if arr[6] == "PASS" and len(alt_lst) == 1:
                alt = alt_lst[0]
                if len(ref) == 1 and len(alt) == 1:
                    hsh[chrom] += 1
    elif vcf_file.endswith(".bgz"):
        tb = tabix.open(vcf_file)
        for chrom in chrom_lst:
            try:
                records = tb.query(chrom, 0, tname2tsize[chrom])
                hsh[chrom] = len(list(records))
            except tabix.TabixError:
                continue

    state = 0
    for chrom in chrom_lst:
        if hsh[chrom] == 0:
            print("{} does not have any variants on chrom: {}".format(vcf_file, chrom))
            state = 1
    if state == 1:
        print("{} might be corrupted".format(vcf_file))
        return 1
    else:
        return 0


def check_vcf_file(
    vcf_file: str, chrom_lst: List[str], tname2tsize: Dict[str, int]
):
    if vcf_file is None:
        print("Please provide the path to the VCF file for the following argument: --vcf")
        return 1
    else:
        if os.path.exists(vcf_file):
            if vcf_file.endswith(".vcf"):
                if os.path.getsize(vcf_file) == 0:
                    return 1
                else:
                    if is_vcf_file_corrupt(vcf_file, chrom_lst, tname2tsize):
                        return 0
                    return 0
            elif vcf_file.endswith(".bgz"):
                tbi_file = vcf_file + ".tbi"
                if os.path.exists(tbi_file):
                    if is_vcf_file_corrupt(vcf_file, chrom_lst, tname2tsize):
                        return 0
                    return 0
                else:
                    print(
                        "tabix index file does not exist for {}".format(vcf_file)
                    )
                    return 1
            elif vcf_file.endswith(".gz"):
                print("hapfusion doesn't support loading of gzip compressed VCF files {}".format(vcf_file))
                return 1
            else:
                print(
                    "VCF file must have the .vcf suffix {}".format(vcf_file)
                )
                return 1
        else:
            print("{} file is missing".format(vcf_file))
            return 1


def is_recomb_file_corrupt(
    recomb_file: str, 
    chrom_lst: List[str], 
):

    hsh = defaultdict(lambda: 0)
    for line in open(recomb_file):
        if line.startswith("#"):
            continue
        fusion = FUSION(line)
        if fusion.is_pass: 
            hsh[fusion.tname] += 1

    state = 0
    for chrom in chrom_lst:
        if hsh[chrom] == 0:
            print("{} does not have any variants on chrom: {}".format(recomb_file, chrom))
            state = 1
            
    if state == 1:
        print("{} might be corrupted".format(recomb_file))
        return 1
    else:
        return 0


def check_recomb_file(
    recomb_file: str, 
    chrom_lst: List[str], 
):
    if recomb_file is None:
        print("Please provide the path to the VCF file for the following argument: --fusion")
        return 1
    else:
        if os.path.exists(recomb_file):
            if os.path.getsize(recomb_file) == 0:
                return 1
            else:
                if is_recomb_file_corrupt(recomb_file, chrom_lst):
                    return 0
                return 0
        else:
            print("{} file is missing".format(recomb_file))
            return 1


def check_out_file(out_file: str):
    if out_file is None:
        return 1
    else:
        if out_file.endswith(".gz"):
            print("hapfusion doesn't support return gzip compressed VCF files")
            return 1
        elif out_file.endswith(".bgz"):
            print("hapfusion doesn't support return bgzip compressed VCF files")
            return 1
        else:
            return 0


def check_phased_vcf_file(
    vcf_file: str, chrom_lst: List[str], tname2tsize: Dict[str, int]
):
    if vcf_file is None:
        print(
            "Please provide the path to the VCF file for the following argument: --phased_vcf"
        )
        return 1
    else:
        if os.path.exists(vcf_file):
            if vcf_file.endswith(".vcf"):
                if os.path.getsize(vcf_file) == 0:
                    return 1
                else:
                    vcf_header_lst = []
                    for line in open(vcf_file).readlines():
                        if line.startswith("#"):
                            vcf_header_lst.append(line.strip())
                        if line.startswith("#CHROM"):
                            break

                    if any(
                        "##FORMAT=<ID=PS" in line.strip() for line in vcf_header_lst
                    ):
                        if is_vcf_file_corrupt(
                            "--phased_vcf", vcf_file, chrom_lst, tname2tsize
                        ):
                            return 1
                        else:
                            return 0
                    else:
                        print("VCF file is not phased")
                        return 1

            elif vcf_file.endswith(".bgz"):
                tbi_file = vcf_file + ".tbi"
                if os.path.exists(tbi_file):
                    v = cyvcf2.VCF(vcf_file)
                    vcf_header_lst = v.raw_header.strip().split()
                    if any(
                        "##FORMAT=<ID=PS" in line.strip() for line in vcf_header_lst
                    ):
                        if is_vcf_file_corrupt(
                            "--phased_vcf", vcf_file, chrom_lst, tname2tsize
                        ):
                            return 1
                        else:
                            return 0
                    else:
                        print("VCF file is not phased")
                        return 1
                else:
                    print("tabix index file does not exist")
                    return 1
            elif vcf_file.endswith(".gz"):
                print("hapfusion doesn't support loading of gzip compressed VCF files")
                return 1
            else:
                print("VCF file must have the .vcf suffix")
                return 1
        else:
            print("{} file is missing".format(vcf_file))
            return 1



def check_caller_input_exists(
    bam_file: str,
    vcf_file: str,
    out_file: str,
    chrom_lst: List[str],
    tname2tsize: Dict[str, int]
) -> None:

    counter = 0
    counter += check_bam_file(bam_file, chrom_lst)
    counter += check_vcf_file(vcf_file, chrom_lst, tname2tsize)
    counter += check_out_file(out_file)
    if counter > 0:
        print("One or more inputs and parameters are missing")
        print("Please provide the correct inputs and parameters")
        exit()


def check_phaser_input_exists(
    bam_file: str,
    vcf_file: str,
    out_file: str,
    chrom_lst: List[str],
    tname2tsize: Dict[str, int],
):

    counter = 0
    counter += check_bam_file(bam_file, chrom_lst)
    counter += check_vcf_file(vcf_file, chrom_lst, tname2tsize)
    counter += check_out_file(out_file)
    if not out_file.endswith(".phased.vcf"):
        print("Please use the suffix .phased.vcf for the output file")
        counter += 1

    if counter > 0:
        print("One or more inputs and parameters are missing")
        print("Please provide the correct inputs and parameters")
        exit()


def check_plot_input_exists(
    bam_file: str,
    vcf_file: str,
    recomb_file: str,
    chrom_lst: List[str],
    tname2tsize: Dict[str, int],
):

    counter = 0
    counter += check_bam_file(bam_file, chrom_lst)
    counter += check_recomb_file(recomb_file, chrom_lst)
    counter += check_vcf_file(vcf_file, chrom_lst, tname2tsize)
    if counter > 0:
        print("One or more inputs and parameters are missing")
        print("Please provide the correct inputs and parameters")
        exit()


