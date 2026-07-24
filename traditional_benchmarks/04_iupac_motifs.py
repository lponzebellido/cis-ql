#!/usr/bin/env python3
"""
04_iupac_motifs.py
Conventional Python script converting IUPAC degenerate ambiguity codes into regular expressions.
Equivalent Cis-QL query: FIND MOTIF "TATAWAW" STRAND POSITIVE AS tata_boxes;
"""

import re

IUPAC_MAP = {
    'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T',
    'R': '[AG]', 'Y': '[CT]', 'S': '[GC]', 'W': '[AT]',
    'K': '[GT]', 'M': '[AC]', 'B': '[CGT]', 'D': '[AGT]',
    'H': '[ACT]', 'V': '[ACG]', 'N': '[ACGT]'
}

def iupac_to_regex(pattern):
    regex = ''
    for char in pattern.upper():
        regex += IUPAC_MAP.get(char, char)
    return regex

def find_iupac_motif(fasta_seq, motif_pattern):
    regex_str = iupac_to_regex(motif_pattern)
    compiled = re.compile(regex_str)
    matches = []
    for match in compiled.finditer(fasta_seq):
        matches.append((match.start(), match.end(), match.group()))
    return matches

if __name__ == '__main__':
    seq = "ATGCGATATATATAGCTAGCTATATATAGCGATCGATCG"
    matches = find_iupac_motif(seq, "TATAWAW")
    print(f"Found {len(matches)} IUPAC motif hits using pattern '{iupac_to_regex('TATAWAW')}'")
