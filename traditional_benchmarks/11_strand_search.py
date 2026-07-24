#!/usr/bin/env python3
"""
11_strand_search.py
Conventional Python script for strand-specific motif searching (sense vs antisense).
Equivalent Cis-QL Query:
FIND MOTIF "TATAAT" STRAND POSITIVE AS sense_promoters;
FIND MOTIF "TATAAT" STRAND NEGATIVE AS antisense_promoters;
EXTRACT sense_promoters;
"""

import re

def reverse_complement(seq):
    trans = str.maketrans('ATGCatgc', 'TACGtacg')
    return seq.translate(trans)[::-1]

def strand_specific_search(sequence, motif="TATAAT"):
    sense_matches = [m.start() for m in re.finditer(motif, sequence)]
    rev_seq = reverse_complement(sequence)
    antisense_matches = [len(sequence) - m.end() for m in re.finditer(motif, rev_seq)]
    return sense_matches, antisense_matches

if __name__ == '__main__':
    seq = "TATAATGCGATCGATATTATA"
    sense, antisense = strand_specific_search(seq)
    print(f"Found {len(sense)} sense and {len(antisense)} antisense motif matches.")
