#!/usr/bin/env python3
"""
08_ctcf_insulators.py
Conventional Python script scanning CTCF PWM and using EXCEPT interval subtraction to remove CDS overlaps.
Equivalent Cis-QL Query: SCAN ctcf_matrix THRESHOLD 75 % AS ctcf_sites; EXCEPT ctcf_sites FROM CDS;
"""

def except_set_subtraction(ctcf_sites, cds_regions):
    intergenic_ctcf = []
    for c_start, c_end in ctcf_sites:
        overlaps = False
        for g_start, g_end in cds_regions:
            if max(c_start, g_start) <= min(c_end, g_end):
                overlaps = True
                break
        if not overlaps:
            intergenic_ctcf.append((c_start, c_end))
    return intergenic_ctcf

if __name__ == '__main__':
    ctcf_sites = [(100, 119), (500, 519), (1000, 1019)]
    cds_regions = [(450, 600)]
    intergenic = except_set_subtraction(ctcf_sites, cds_regions)
    print(f"Isolated {len(intergenic)} intergenic CTCF sites after EXCEPT subtraction.")
