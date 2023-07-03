import math
import pysam
import random
import natsort
import numpy as np
import himut.cslib
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


    def cs2tpos2qbase(self):

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

    def get_tcoord(self):
        self.tcoord = "{}:{}-{}".format(self.tname, self.tstart, self.tend)

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

