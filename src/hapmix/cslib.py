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


def cs2tpos2qbase(read) -> Dict[int, Tuple[str, int]]:
    """
    Converts cstuple to 1-coordinate based read allele
    
    Parameters:
        tpos (int): reference alignemnt start position
        qpos (int): query alignment start position
        qbq_lst: list containing base quality scores for CCS bases

    Returns:
        dictionary mapping reference position to tuple containing reference base, alternative base and base quality score
    """

    tpos = read.tstart
    qpos = read.qstart
    read.tpos2qpos = {}
    read.tpos2qbase = {}
    for cstuple in read.cstuple_lst:
        state, ref, alt, ref_len, alt_len, = cstuple
        if state == 1:  # match
            for i, alt_base in enumerate(alt):
                read.tpos2qpos[tpos + i + 1] = qpos + i
                read.tpos2qbase[tpos + i + 1] = (alt_base, read.bq_int_lst[qpos + i])
        elif state == 2:  # sub
            read.tpos2qpos[tpos + 1] = qpos
            read.tpos2qbase[tpos + 1] = (alt, read.bq_int_lst[qpos])
        elif state == 3:  # insertion ## to do for hetINDEL phasing
            pass
        elif state == 4:  # deletion ## to do for hetINDEL phasing
            for j in range(len(ref[1:])):
                read.tpos2qbase[tpos + j + 1] = ("-", 0)
        tpos += ref_len
        qpos += alt_len


def cs2mismatch(read): ## todo

    state = 0
    tpos = read.tstart
    mismatch_lst = []
    for cstuple in read.cstuple_lst:
        mstate, ref, alt, ref_len, alt_len, = cstuple
        if state == 0 and mstate == 1:  # init # match
            state = 0
            counter = 0
        elif state == 0 and mstate != 1:  # init # mismatch
            counter += 1
            ref_lst = [ref]
            alt_lst = [alt]
            if mstate == 2:  # snp
                state = 1
                tstart = tpos
            elif mstate == 3 or mstate == 4:  # insertion # deletion
                state = 2
                tstart = tpos - 1
        elif state != 0 and mstate == 2:  # snp
            state = 1
            counter += 1
            ref_lst.append(ref)
            alt_lst.append(alt)
        elif state != 0 and mstate == 3:  # insertion
            state = 2
            counter += 1
            ref_lst.append(ref)
            alt_lst.append(alt)
        elif state != 0 and mstate == 4:  # deletion
            state = 2
            counter += 1
            ref_lst.append(ref)
            alt_lst.append(alt)
        elif (
            state != 0 and mstate == 1 and ref_len <= 10
        ):  # match # mnp: condition # snp, match, snp
            counter += 1
            state = state
            ref_lst.append(ref)
            alt_lst.append(ref)
        elif state != 0 and mstate == 1 and ref_len > 11:  # match # return
            state = 4
        tpos += ref_len 

        # return
        if state == 4:
            ref = "".join(ref_lst)
            alt = "".join(alt_lst)
            if len(ref) == 1 and len(alt) == 1: # snv
                if ref == "N":
                    continue
                mismatch_lst.append((read.tname, tstart + 1, ref, alt))
            else:
                mismatch_lst.append((read.tname, tstart + 1, ref, alt))
            state = 0  
    read.mismatch_lst = mismatch_lst
