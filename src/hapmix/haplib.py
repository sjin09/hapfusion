import bisect
from collections import defaultdict
from typing import Dict, List, Tuple


def get_hamming_distance(consensus_h0_lst: List[str], read_hbit_lst: List[str]) -> Tuple[int, int]:

    h0_distance = 0
    h1_distance = 0
    for (i, j) in zip(consensus_h0_lst, read_hbit_lst):
        if j == "-":
            continue
        if i != j:
            h0_distance += 1
        else:
            h1_distance += 1
    return h0_distance, h1_distance


def get_read_haplotype(consensus_h0_lst: List[str], read_hbit_lst: List[str]) -> str: 
    h0_distance, h1_distance = get_hamming_distance(consensus_h0_lst, read_hbit_lst)
    if h0_distance > h1_distance:
        return "1"
    elif h0_distance < h1_distance:
        return "0"
    elif h0_distance == h1_distance:
        return "."


def get_hapmix_idx(consensus_h0_lst, read_hbit_lst):
    idx_lst = []
    read_haplotype = get_read_haplotype(consensus_h0_lst, read_hbit_lst)
    if read_haplotype == "0":
        for i, (j, k) in enumerate(zip(consensus_h0_lst, read_hbit_lst)):
            if k == "-":
                continue
            if j != k:
                idx_lst.append(i)
    elif read_haplotype == "1":
        for x, (y, z) in enumerate(zip(consensus_h0_lst, read_hbit_lst)):
            if z == "-":
                continue
            if y == z:
                idx_lst.append(x)
    elif read_haplotype == ".":
        for i, (j, k) in enumerate(zip(consensus_h0_lst, read_hbit_lst)):
            if k == "-":
                continue
            if j != k:
                idx_lst.append(i)
    return idx_lst, read_haplotype



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


def get_read_hetsnps(
    tpos2qbase: Dict[int, Tuple[str, int]],
    hpos_lst: List[int],  
    hetsnp_lst: List[Tuple[str, int, str, str]],
) -> List[Tuple[str, int, str, str]]:

    tpos_lst = list(tpos2qbase.keys())
    idx = bisect.bisect_right(hpos_lst, tpos_lst[0])
    jdx = bisect.bisect_right(hpos_lst, tpos_lst[-1])
    if idx == jdx:
        return 0, []
    elif jdx - idx == 1:
        return 1, []
    else:
        phased_hetsnp_subset_lst = [hetsnp_lst[kdx] for kdx in range(idx, jdx)]
    return len(phased_hetsnp_subset_lst), phased_hetsnp_subset_lst


def get_read_hbit_lst(
    tpos2qbase: Dict[int, Tuple[str, int]],
    hetsnp2bidx: Dict[Tuple[str, int, str, str], int],
    hetsnp2hstate: Dict[Tuple[str, int, str, str], str],
    hetsnp_subset_lst: List[Tuple[str, int, str, str]],
):

    bidx2hetsnp_subset_lst = defaultdict(list)
    for hetsnp in hetsnp_subset_lst:
        bidx = hetsnp2bidx[hetsnp]
        bidx2hetsnp_subset_lst[bidx].append(hetsnp)
    cnt = len(bidx2hetsnp_subset_lst.keys()) ## todo
    
    if cnt == 1:
        hetsnp_subset_lst = 0 
        hetsnp_subset_count = 0
        for _, hetsnp_subset_candidate_lst in bidx2hetsnp_subset_lst.items():
            hetsnp_subset_candidate_count = len(hetsnp_subset_candidate_lst)
            if hetsnp_subset_count < hetsnp_subset_candidate_count:
                hetsnp_subset_lst = hetsnp_subset_candidate_lst
                hetsnp_subset_count = hetsnp_subset_candidate_count

        read_hbit_lst = []
        for hetsnp in hetsnp_subset_lst: # get read haplotype bits
            qbase, _ = tpos2qbase[hetsnp[1]]
            if qbase == hetsnp[2]:
                read_hbit_lst.append("0")
            elif qbase == hetsnp[3]:
                read_hbit_lst.append("1")
            else:
                read_hbit_lst.append("-")

        h0_hbit_lst = [hetsnp2hstate[hetsnp] for hetsnp in hetsnp_subset_lst]
        h1_hbit_lst = get_haplotype_complement(h0_hbit_lst)
        return read_hbit_lst, h0_hbit_lst, h1_hbit_lst 
    else:
        return [], [], [] 
