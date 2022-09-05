import bisect
from collections import defaultdict
from typing import Dict, List, Tuple


def get_hamming_distance(
    h0_hbit_lst: List[str], 
    ccs_hbit_lst: List[str]
) -> Tuple[int, int]:

    h0_distance = 0
    h1_distance = 0
    for (i, j) in zip(h0_hbit_lst, ccs_hbit_lst):
        if j == "-":
            continue
        if i != j:
            h0_distance += 1
        else:
            h1_distance += 1
    return h0_distance, h1_distance


def get_ccs_haplotype(
    h0_hbit_lst: List[str], 
    ccs_hbit_lst: List[str]
) -> str: 

    h0_distance, h1_distance = get_hamming_distance(h0_hbit_lst, ccs_hbit_lst)
    if h0_distance > h1_distance:
        return "1"
    elif h0_distance < h1_distance:
        return "0"
    elif h0_distance == h1_distance:
        return "."


def get_haplotype_complement(hbit_lst: List[str]) -> List[str]:
    hbit_complement_lst = []
    for hbit in hbit_lst:
        if hbit == "0":
            hbit_complement_lst.append("1")
        elif hbit == "1":
            hbit_complement_lst.append("0")
        else:
            hbit_complement_lst.append(hbit)
    return hbit_complement_lst


def get_ccs_hbit_lst(
    ccs,
    hetsnp2bidx: Dict[Tuple[str, int, str, str], int],
    hetsnp2hstate: Dict[Tuple[str, int, str, str], str],
    phased_hpos_lst: List[Tuple[str, int, str, str]],
    phased_hetsnp_lst: List[Tuple[str, int, str, str]],
):

    idx = bisect.bisect_right(phased_hpos_lst, ccs.tstart)
    jdx = bisect.bisect_right(phased_hpos_lst, ccs.tend)
    if idx == jdx:
        return 0, []
    elif jdx - idx == 1:
        return 1, []
    else:
        phased_hetsnp_subset_lst = [phased_hetsnp_lst[kdx] for kdx in range(idx, jdx)]

    bidx2hetsnp_subset_lst = defaultdict(list)
    for phased_hetsnp in phased_hetsnp_subset_lst:
        bidx = hetsnp2bidx[phased_hetsnp]
        bidx2hetsnp_subset_lst[bidx].append(phased_hetsnp)
    blk_cnt = len(bidx2hetsnp_subset_lst.keys()) ## todo
   
    if blk_cnt == 1:
        h0_hbit_lst = []
        ccs_hbit_lst = []
        blk_hetsnp_lst = list(bidx2hetsnp_subset_lst.values())[0]
        for hetsnp in blk_hetsnp_lst: # get read haplotype bits
            qbase, _ = ccs.tpos2qbase[hetsnp[0]]
            h0_hbit_lst.append(hetsnp2hstate[hetsnp])
            if qbase == hetsnp[1]: # ref
                ccs_hbit_lst.append("0")
            elif qbase == hetsnp[2]: # alt
                ccs_hbit_lst.append("1")
            else: # del
                ccs_hbit_lst.append("-")
        h1_hbit_lst = get_haplotype_complement(h0_hbit_lst)
        return ccs_hbit_lst, h0_hbit_lst, h1_hbit_lst, phased_hetsnp_subset_lst 
    else:
        return [], [], [], [] 


def get_hapsmash_hetsnps(
    ccs_hap: str,
    h0_hbit_lst: List[str], 
    ccs_hbit_lst: List[str],
    phased_hetsnp_subset_lst: List[Tuple[int, str, str]]
) -> Tuple[List[int], str]:

    hapsmash_hetsnp_lst = []
    hapsmash_hetsnp_idx_lst = []
    if ccs_hap == "0":
        for i, (j, k) in enumerate(zip(h0_hbit_lst, ccs_hbit_lst)):
            if k == "-":
                continue
            if j != k:
                hapsmash_hetsnp_idx_lst.append(i)
                hapsmash_hetsnp_lst.append(phased_hetsnp_subset_lst[i])
    elif ccs_hap == "1":
        for x, (y, z) in enumerate(zip(h0_hbit_lst, ccs_hbit_lst)):
            if z == "-":
                continue
            if y == z:
                hapsmash_hetsnp_idx_lst.append(x)
                hapsmash_hetsnp_lst.append(phased_hetsnp_subset_lst[x])
    else:
        h0_distance = 0
        h1_distance = 0
        for (i, j) in zip(h0_hbit_lst, ccs_hbit_lst):
            if i != j:
                h0_distance += 1
            else:
                h1_distance += 1
        if h0_distance < h1_distance: # ccs haplotype: 0
            for i, (j, k) in enumerate(zip(h0_hbit_lst, ccs_hbit_lst)):
                if k == "-":
                    continue
                if j != k:
                    hapsmash_hetsnp_idx_lst.append(i)
                    hapsmash_hetsnp_lst.append(phased_hetsnp_subset_lst[i])
        else: # ccs haplotype: 1 
            for i, (j, k) in enumerate(zip(h0_hbit_lst, ccs_hbit_lst)):
                if k == "-":
                    continue
                if j == k:
                    hapsmash_hetsnp_idx_lst.append(i)
                    hapsmash_hetsnp_lst.append(phased_hetsnp_subset_lst[i])
    return hapsmash_hetsnp_lst, hapsmash_hetsnp_idx_lst
