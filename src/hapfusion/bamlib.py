import math
import pysam
import bisect
import random
import natsort
import numpy as np
import himut.cslib
import hapfusion.haplib
from typing import Set, List, Dict, Tuple


class BAM:
    def __init__(self, line):
        # target
        if line.is_secondary:
            self.is_primary = False
        else:
            self.is_primary = True
            self.tname = line.reference_name
            self.tstart = line.reference_start
            self.tend = line.reference_end
            # query
            self.qname = line.query_name
            self.qstart = line.query_alignment_start
            self.qend = line.query_alignment_end
            self.qseq = line.query_sequence
            self.qlen = len(self.qseq)
            self.mapq = line.mapping_quality
            self.bq_int_lst = line.query_qualities
            himut.cslib.cs2tuple(self, line.get_tag("cs"))

    def get_qv(self):
        self.qv = np.mean(self.bq_int_lst)
        return self.qv

    def get_tcoord(self):
        self.tcoord = "{}:{}-{}".format(self.tname, self.tstart, self.tend)

    # def get_ccs_hap(
    #     self,
    #     hbit_lst: Dict[str, List[str]],
    #     hpos_lst: Dict[str, List[int]],
    #     hetsnp_lst: List[Tuple[int, str, str]],
    # ):

    #     idx = bisect.bisect_right(hpos_lst, self.tstart)
    #     jdx = bisect.bisect_right(hpos_lst, self.tend)
    #     if (jdx - idx) < 3:
    #         self.hap = "."
    #     else:
    #         self.hetsnp_lst = [hetsnp_lst[kdx] for kdx in range(idx, jdx)]
    #         self.h0_hbit = "".join([hbit_lst[kdx] for kdx in range(idx, jdx)])
    #         self.h1_hbit = "".join([hapfusion.haplib.bit_complement_hsh[hbit] for hbit in self.h0_hbit])
    #         self.hbit = hapfusion.haplib.get_self_hbit(self, self.hetsnp_lst)
    #         if self.h0_hbit == self.hbit:
    #             self.hap = "0"
    #         elif self.h1_hbit == self.hbit:
    #             self.hap = "1"
    #         else:
    #             h0_hd, h1_hd = hapfusion.haplib.get_hap_hamming_distance(self.h0_hbit, self.hbit)
    #             if h0_hd > h1_hd:
    #                 self.hap = "1"
    #             elif h0_hd < h1_hd:
    #                 self.hap = "0"
    #             elif h0_hd == h1_hd:
    #                 self.hap = "2"

    def get_cs2tpos2qbase(self):

        tpos = self.tstart
        qpos = self.qstart
        self.rpos2qpos = {}
        self.tpos2qbase = {}
        self.mismatch_lst = [] 
        for (state, ref, alt, ref_len, alt_len) in self.cstuple_lst:
            if state == 1:  # match
                for i, alt_base in enumerate(alt):
                    self.rpos2qpos[tpos+i] = qpos + i
                    self.tpos2qbase[tpos + i + 1] = (alt_base, self.bq_int_lst[qpos + i])
            elif state == 2:  # sub
                self.rpos2qpos[tpos] = qpos
                self.tpos2qbase[tpos + 1] = (alt, self.bq_int_lst[qpos])
            elif state == 3:  # insertion
                self.mismatch_lst.append((tpos+1, ref, alt))
            elif state == 4:  # deletion
                self.rpos2qpos[tpos] = qpos
                for j in range(len(ref[1:])):
                    self.tpos2qbase[tpos + j + 1] = ("-", 0)
                self.mismatch_lst.append((tpos+1, ref, alt))
            tpos += ref_len
            qpos += alt_len

    def get_blast_sequence_identity(self):
 
        match_count = 0
        mismatch_count = 0
        for cstuple in self.cstuple_lst:
            mstate, _, _, ref_len, alt_len = cstuple
            if mstate == 1:  # match
                match_count += ref_len
            elif mstate == 2:  # mismatch: snp
                mismatch_count += alt_len
            elif mstate == 3:  # mismatch: insertion
                mismatch_count += alt_len
            elif mstate == 4:  # mismatch: deletion
                mismatch_count += ref_len
        alignment_len = match_count + mismatch_count
        blast_sequence_identity = (match_count)/float(alignment_len)
        return blast_sequence_identity

    def get_query_alignment_proportion(self):
        qaln_proportion = (self.qend - self.qstart)/float(self.qlen)
        return qaln_proportion

    def load_mutations(
        self, 
        het_set: Set[Tuple[str, int, str, str]],
        hom_set: Set[Tuple[str, int, str, str]],
        hetpos_set: Set[int],
        hompos_set: Set[int],
        phased_hetsnp_set: Set[Tuple[str, int, str, str]],
        pos2allele: Dict[int, List[str]],
    ):

        self.hetsnp_lst = []
        self.homsnp_lst = []
        self.hetindel_lst = []
        self.homindel_lst = []
        self.denovo_sbs_lst = []
        self.denovo_indel_lst = []
        for (tpos, ccs_ref, ccs_alt, bq) in self.tsbs_lst:
            alpha = bq/float(93)
            # print(tpos, ccs_ref, ccs_alt, bq)
            if (tpos, ccs_ref, ccs_alt) in hom_set:
                self.homsnp_lst.append((tpos, ccs_ref, ccs_alt, alpha))
            elif (tpos, ccs_ref, ccs_alt) in het_set:
                continue
            elif (tpos, ccs_ref, ccs_alt) in phased_hetsnp_set:
                self.hetsnp_lst.append((tpos, ccs_ref, ccs_alt, alpha))
            else:
                self.denovo_sbs_lst.append((tpos, ccs_ref, ccs_alt, alpha))
               
        for (tpos, ccs_ref, ccs_alt, bq) in self.tindel_lst: ## TODO, indels can be approximately similar
            alpha = bq/float((93 * len(ccs_alt)))
            if (tpos, ccs_ref, ccs_alt) in het_set:
                self.hetindel_lst.append((tpos, ccs_ref, ccs_alt, alpha))
            elif (tpos, ccs_ref, ccs_alt) in hom_set:
                self.homindel_lst.append((tpos, ccs_ref, ccs_alt, alpha))
            else:
                if tpos in pos2allele:
                    ccs_ref_len = len(ccs_ref) 
                    ccs_alt_len = len(ccs_alt) 
                    if ccs_ref_len > ccs_alt_len:
                        ccs_state = 0
                    elif ccs_ref_len < ccs_alt_len:
                        ccs_state = 1

                    vcf_ref_len, vcf_alt_len = [len(allele) for allele  in pos2allele[tpos]]
                    if vcf_ref_len > vcf_alt_len:
                        vcf_state = 0
                    elif vcf_ref_len < vcf_alt_len: 
                        vcf_state = 1
                    
                    if tpos in hetpos_set: 
                        if ccs_state == 0 and vcf_state == 0:
                            self.hetindel_lst.append((tpos, ccs_ref, ccs_alt, alpha))
                        elif ccs_state == 1 and vcf_state == 1:
                            self.hetindel_lst.append((tpos, ccs_ref, ccs_alt, alpha))
                        else:
                            self.denovo_indel_lst.append((tpos, ccs_ref, ccs_alt, alpha))
                    elif tpos in hompos_set:
                        if ccs_state == 0 and vcf_state == 0:
                            self.homindel_lst.append((tpos, ccs_ref, ccs_alt, alpha))
                        elif ccs_state == 1 and vcf_state == 1:
                            self.homindel_lst.append((tpos, ccs_ref, ccs_alt, alpha))
                        else:
                            self.denovo_indel_lst.append((tpos, ccs_ref, ccs_alt, alpha))
                else:  
                    self.denovo_indel_lst.append((tpos, ccs_ref, ccs_alt, alpha))
                    
        self.mismatch_lst = natsort.natsorted(self.denovo_sbs_lst + self.denovo_indel_lst + self.hetindel_lst + self.homindel_lst)
        self.mismatch_tpos_lst = [mismatch[0] for mismatch in self.mismatch_lst]
 

def get_sample(bam_file: str):

    state = 0
    alignments = pysam.AlignmentFile(bam_file, "rb")
    bam_header_lines = str(alignments.header).strip().split("\n")
    for line in bam_header_lines:
        if line.startswith("@RG"):
            for i in line.split():
                if i.startswith("SM"):
                    sample = i.split(":")[1]
                    return sample
    if state == 0:
        print("SM field is missing")
        print("Please provide a BAM file with @RG group")
        print(
            "samtools reheader in.header.sam in.bam > out.bam command can be used to insert a new header"
        )
        hapfusion.util.exit()


def get_tname2tsize(bam_file: str) -> Tuple[List[str], Dict[str, int]]:

    tname2tsize = {}
    alignments = pysam.AlignmentFile(bam_file, "rb")
    bam_header_lst = str(alignments.header).strip().split("\n")
    for h in bam_header_lst:
        if h.startswith("@SQ"):
            arr = h.split("\t")
            tname = arr[1].replace("SN:", "")
            tsize = int(arr[2].replace("LN:", ""))
            tname2tsize[tname] = tsize
    alignments.close()

    tname_lst = natsort.natsorted(list(tname2tsize.keys()))
    if len(tname_lst) == 0:
        print("@SQ header is missing from BAM file")
        print(
            "Please use samtools reheader to insert approprirate header to your BAM file"
        )
        himut.util.exit()
    return tname_lst, tname2tsize


def get_md_threshold(coverage: int) -> int:
    md_threshold = math.ceil(coverage + (4 * math.sqrt(coverage)))
    return md_threshold


def get_thresholds(
    bam_file: str,
    chrom_lst: List[str],
    chrom2len: Dict[str, int],
) -> Tuple[int, int, int, int]:

    if len(chrom_lst) == 0:
        print("target is missing")
        print("Please check .vcf file or .target file")
        himut.util.exit()

    qlen_lst = []
    random.seed(10)
    sample_count = 100
    sample_range = 100000
    genome_read_sum = 0
    genome_sample_sum = sample_count * sample_range * len(chrom_lst)
    alignments = pysam.AlignmentFile(bam_file, "rb")
    for chrom in chrom_lst:
        chrom_len = chrom2len[chrom]
        random_start_lst = random.sample(range(chrom_len), sample_count)
        for start in random_start_lst:
            end = start + 100000
            for read in alignments.fetch(chrom, start, end):
                mapq = int(read.mapping_quality)
                alignment_type = read.get_tag("tp") if read.has_tag("tp") else "."
                if mapq > 0 and alignment_type == "P":
                    qlen = len(read.query_sequence)
                    genome_read_sum += qlen
                    qlen_lst.append(qlen)
    alignments.close()
    qlen_std = np.std(qlen_lst)
    qlen_mean = math.ceil(np.mean(qlen_lst))
    qlen_lower_limit = (
        0
        if math.ceil(qlen_mean - 2 * qlen_std) < 0
        else math.ceil(qlen_mean - 2 * qlen_std)
    )
    qlen_upper_limit = math.ceil(qlen_mean + 2 * qlen_std)
    coverage = genome_read_sum / float(genome_sample_sum)
    md_threshold = get_md_threshold(coverage)
    return qlen_lower_limit, qlen_upper_limit, md_threshold 
 
  
def get_read_depth(
    basecounts: Dict[int, Dict[int, int]],
):
    ins_count = basecounts[4]
    del_count = basecounts[5]
    indel_count = ins_count + del_count
    read_depth = sum(basecounts) - ins_count
    return read_depth, indel_count


def get_trim_range(
    qlen: int,
    min_trim: float,
):
    trim_qstart = math.floor(min_trim * qlen)
    trim_qend = math.ceil((1 - min_trim) * qlen)
    return trim_qstart, trim_qend
     

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


def get_mismatch_positions(ccs):
    mismatch_tpos_lst = [mismatch[0] for mismatch in ccs.mismatch_lst] # 1-coordinate
    return mismatch_tpos_lst 

