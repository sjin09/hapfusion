import re
from typing import Dict, List, Tuple


def cs2lst(cs_tag):
    cslst = [cs for cs in re.split("(:[0-9]+|\*[a-z][a-z]|[=\+\-][A-Za-z]+)", cs_tag)]
    cslst = [cs.upper() for cs in cslst if cs != ""]
    return cslst


def cs2tuple(read) -> List[Tuple[int, str, str, int, int]]:

    qpos = read.qstart
    read.cstuple_lst = []
    cs_lst = cs2lst(read.cs_tag)
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
            m = read.qseq[qstart:qend]
            t = (1, m, m, mlen, mlen)
        elif cs.startswith("*"):  # snp # target and query
            mlen = 1
            ref, alt = list(m)
            t = (2, ref, alt, 1, 1)
        elif cs.startswith("+"):  # insertion # query
            ref = read.qseq[qpos - 1]
            alt = ref + m
            t = (3, ref, alt, 0, mlen)
        elif cs.startswith("-"):  # deletion # target
            alt = read.qseq[qpos - 1]
            ref = alt + m
            t = (4, ref, alt, mlen, 0)
            mlen = 0
        read.cstuple_lst.append(t)
        qpos += mlen


def cs2tpos2qbase(ccs) -> Dict[int, Tuple[str, int]]:
    """
    Converts cstuple to 1-coordinate based read allele
    
    Parameters:
        tpos (int): reference alignemnt start position
        qpos (int): query alignment start position
        qbq_lst: list containing base quality scores for CCS bases

    Returns:
        dictionary mapping reference position to tuple containing reference base, alternative base and base quality score
    """

    tpos = ccs.tstart
    qpos = ccs.qstart
    ccs.tpos2qpos = {}
    ccs.tpos2qbase = {}
    for cstuple in ccs.cstuple_lst:
        state, ref, alt, ref_len, alt_len, = cstuple
        if state == 1:  # match
            for i, alt_base in enumerate(alt):
                ccs.tpos2qpos[tpos + i + 1] = qpos + i
                ccs.tpos2qbase[tpos + i + 1] = (alt_base, ccs.bq_int_lst[qpos + i])
        elif state == 2:  # sub
            ccs.tpos2qpos[tpos + 1] = qpos
            ccs.tpos2qbase[tpos + 1] = (alt, ccs.bq_int_lst[qpos])
        elif state == 3:  # insertion ## to do for hetINDEL phasing
            pass
        elif state == 4:  # deletion ## to do for hetINDEL phasing
            for j in range(len(ref[1:])):
                ccs.tpos2qbase[tpos + j + 1] = ("-", 0)
        tpos += ref_len
        qpos += alt_len


def cs2subindel(ccs):

    tpos = ccs.tstart
    qpos = ccs.qstart
    ccs.tsbs_lst = []
    ccs.qsbs_lst = []
    ccs.qsbs_bq_lst = []
    ccs.tdel_lst = []
    ccs.tins_lst = []
    ccs.qins_bq_lst = []
    for cstuple in ccs.cstuple_lst:
        state, ref, alt, ref_len, alt_len, = cstuple
        if ref.count("N") > 0:
            continue
        if state == 2:  # snp 
            ccs.qsbs_bq_lst.append(ccs.bq_int_lst[qpos])
            ccs.tsbs_lst.append((ccs.tname, tpos + 1, ref, alt))
            ccs.qsbs_lst.append((ccs.qname, qpos + 1, alt, ref))
        elif state == 3: # insertion
            ins_start = qpos - 1
            ins_end = qpos - 1 + len(alt)
            ccs.tins_lst.append((ccs.tname, tpos, ref, alt))
            ccs.qins_bq_lst.append(ccs.bq_int_lst[ins_start:ins_end])
        elif state == 4: # deletion
            ccs.tdel_lst.append((ccs.tname, tpos, ref, alt))
        tpos += ref_len 
        qpos += alt_len


