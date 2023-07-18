import bisect
from typing import Dict, List, Tuple
hap_lst = ["0", "1", "2"] 
bit_complement_hsh = {"0": "1", "1": "0", "-": "-"}


def get_ccs_hbit(ccs, hetsnp_subset_lst: List[Tuple[str, int, str, str]]) -> List[str]:

    ccs_hbit = ""
    ccs_hbit_bq_lst = []
    for hetsnp in hetsnp_subset_lst:  # get read haplotype bits
        qbase, qbq = ccs.tpos2qbase[hetsnp[0]]
        if qbase == hetsnp[1]:  # ref
            ccs_hbit += "0"
            ccs_hbit_bq_lst.append(qbq)
        elif qbase == hetsnp[2]:  # alt
            ccs_hbit += "1"
            ccs_hbit_bq_lst.append(qbq)
        else:
            if qbase == "-":
                ccs_hbit += "-"
                ccs_hbit_bq_lst.append("!")
            else:
                ccs_hbit += "2"
                ccs_hbit_bq_lst.append(qbq)
            
    return ccs_hbit


def get_hap_hamming_distance(
    h0_hbit: str, 
    hbit: str
):

    h0_distance = 0
    h1_distance = 0
    for (i, j) in zip(h0_hbit, hbit):
        if j == "-":
            continue
        elif j == "2":
            continue
        if i != j:
            h0_distance += 1
        else:
            h1_distance += 1
    return h0_distance, h1_distance


def get_ccs_hap(
    ccs,
    hbit_lst: Dict[str, List[str]],
    hpos_lst: Dict[str, List[int]],
    hetsnp_lst: List[Tuple[int, str, str]],
):

    idx = bisect.bisect_right(hpos_lst, ccs.tstart)
    jdx = bisect.bisect_right(hpos_lst, ccs.tend)
    if (jdx - idx) < 3:
        ccs.hap = "."
    else:
        ccs.hetsnp_lst = [hetsnp_lst[kdx] for kdx in range(idx, jdx)]
        ccs.h0_hbit = "".join([hbit_lst[kdx] for kdx in range(idx, jdx)])
        ccs.h1_hbit = "".join([bit_complement_hsh[hbit] for hbit in ccs.h0_hbit])
        ccs.hbit = get_ccs_hbit(ccs, ccs.hetsnp_lst)
        if ccs.h0_hbit == ccs.hbit:
            ccs.hap = "0"
        elif ccs.h1_hbit == ccs.hbit:
            ccs.hap = "1"
        else:
            h0_hd, h1_hd = get_hap_hamming_distance(ccs.h0_hbit, ccs.hbit)
            if h0_hd > h1_hd:
                ccs.hap = "1'"
            elif h0_hd < h1_hd:
                ccs.hap = "0'"
            elif h0_hd == h1_hd:
                ccs.hap = "2"


# def get_hapfusion_hetsnps(
#     ccs_hap: str,
#     h0_hbit_lst: List[str], 
#     ccs_hbit_lst: List[str],
#     phased_hetsnp_subset_lst: List[Tuple[int, str, str]]
# ) -> Tuple[List[int], str]:

#     haphunter_hetsnp_lst = []
#     haphunter_indel_idx_set = set()
#     haphunter_hetsnp_idx_set = set()
#     if ccs_hap == "0":
#         for i, (j, k) in enumerate(zip(h0_hbit_lst, ccs_hbit_lst)):
#             if k == "-":
#                 haphunter_indel_idx_set.add(i)
#                 continue
#             if j != k:
#                 haphunter_hetsnp_idx_set.add(i)
#                 haphunter_hetsnp_lst.append(phased_hetsnp_subset_lst[i])
#     elif ccs_hap == "1":
#         for i, (j, k) in enumerate(zip(h0_hbit_lst, ccs_hbit_lst)):
#             if k == "-":
#                 haphunter_indel_idx_set.add(i)
#                 continue
#             if j == k:
#                 haphunter_hetsnp_idx_set.add(i)
#                 haphunter_hetsnp_lst.append(phased_hetsnp_subset_lst[i])
#     else:
#         h0_distance = 0
#         h1_distance = 0
#         for (i, j) in zip(h0_hbit_lst, ccs_hbit_lst):
#             if i != j:
#                 h0_distance += 1
#             else:
#                 h1_distance += 1
#         if h0_distance < h1_distance: # ccs haplotype: 0 and "."
#             for i, (j, k) in enumerate(zip(h0_hbit_lst, ccs_hbit_lst)):
#                 if k == "-":
#                     haphunter_indel_idx_set.add(i)
#                     continue
#                 if j != k:
#                     haphunter_hetsnp_idx_set.add(i)
#                     haphunter_hetsnp_lst.append(phased_hetsnp_subset_lst[i])
#         else: # ccs haplotype: 1 
#             for i, (j, k) in enumerate(zip(h0_hbit_lst, ccs_hbit_lst)):
#                 if k == "-":
#                     haphunter_indel_idx_set.add(i)
#                     continue
#                 if j == k:
#                     haphunter_hetsnp_idx_set.add(i)
#                     haphunter_hetsnp_lst.append(phased_hetsnp_subset_lst[i])
#     return haphunter_hetsnp_lst, haphunter_hetsnp_idx_set, haphunter_indel_idx_set 
