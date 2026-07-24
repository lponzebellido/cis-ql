#!/usr/bin/env python3
"""
13_complex_where.py
Conventional Python script combining physical length filters and Smith-Waterman similarity scoring.
Equivalent Cis-QL Query: EXTRACT GENE WHERE LENGTH > 800 BP AND SIMILARITY > 65 %;
"""

def smith_waterman_score(seq1, seq2):
    m, n = len(seq1), len(seq2)
    dp = [[0] * (n + 1) for _ in range(m + 1)]
    max_score = 0
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            s = 2 if seq1[i-1] == seq2[j-1] else -1
            score = max(0, dp[i-1][j-1] + s, dp[i-1][j] - 2, dp[i][j-1] - 2)
            dp[i][j] = score
            if score > max_score: max_score = score
    max_possible = min(m, n) * 2
    return (max_score / max_possible * 100.0) if max_possible > 0 else 0.0

def complex_where_filter(genes, ref_seq, min_len=800, min_sim=65.0):
    filtered = []
    for g in genes:
        g_len = g['end'] - g['start'] + 1
        if g_len > min_len:
            sim = smith_waterman_score(g['sequence'], ref_seq)
            if sim > min_sim:
                filtered.append((g, sim))
    return filtered

if __name__ == '__main__':
    print("Complex multi-property filtering script executed.")
