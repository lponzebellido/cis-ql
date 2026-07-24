#!/usr/bin/env python3
"""
10_promoter_union.py
Conventional Python script finding -35 and -10 boxes and executing interval UNION merging.
Equivalent Cis-QL Query:
FIND MOTIF "TTGACA" WITHIN 200 BP UPSTREAM FROM CDS AS minus35_box;
FIND MOTIF "TATAAT" WITHIN 50 BP UPSTREAM FROM CDS AS minus10_box;
UNION minus35_box AND minus10_box;
"""

def union_interval_sets(set1, set2):
    combined = sorted(set1 + set2, key=lambda x: x[0])
    if not combined: return []
    merged = [combined[0]]
    for current in combined[1:]:
        prev_start, prev_end = merged[-1]
        curr_start, curr_end = current
        if curr_start <= prev_end:
            merged[-1] = (prev_start, max(prev_end, curr_end))
        else:
            merged.append(current)
    return merged

if __name__ == '__main__':
    m35 = [(100, 106), (300, 306)]
    m10 = [(104, 110), (500, 506)]
    union_res = union_interval_sets(m35, m10)
    print(f"Merged {len(m35) + len(m10)} sites into {len(union_res)} unified promoter intervals.")
