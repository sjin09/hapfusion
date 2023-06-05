import math
import pysam
import random
import natsort
import numpy as np
import himut.cslib
import hapfusion.util
import hapfusion.vcflib
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

    def cs2mut(self):
        himut.cslib.cs2mut(self)

    def cs2subindel(self):
        himut.cslib.cs2subindel(self)

    def cs2tpos2qbase(self):
        himut.cslib.cs2tpos2qbase(self)

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
   
        