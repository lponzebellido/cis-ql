#!/usr/bin/env python3
"""
03_similarity_align.py
Conventional Python implementation using Smith-Waterman local alignment for similarity filtering.
Equivalent Cis-QL query: EXTRACT GENE WHERE LENGTH > 1 KB AND SIMILARITY > 70 %;
"""

import sys

def smith_waterman(seq1, seq2, match=2, mismatch=-1, gap=-2):
    m, n = len(seq1), len(seq2)
    score_matrix = [[0] * (n + 1) for _ in range(m + 1)]
    max_score = 0
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            s = match if seq1[i-1] == seq2[j-1] else mismatch
            score = max(0, score_matrix[i-1][j-1] + s, score_matrix[i-1][j] + gap, score_matrix[i][j-1] + gap)
            score_matrix[i][j] = score
            if score > max_score:
                max_score = score
    max_possible = min(len(seq1), len(seq2)) * match
    return (max_score / max_possible * 100.0) if max_possible > 0 else 0.0

if __name__ == '__main__':
    ref = "ATGCGATCGATCGATCGATCGATCGATCGATCGATCGATC"
    query = "ATGCGATCGATCGATCGATCGATCGATC"
    sim = smith_waterman(ref, query)
    print(f"Calculated similarity: {sim:.2f}%")
