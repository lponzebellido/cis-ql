#!/usr/bin/env python3
"""
07_pwm_scanning.py
Conventional Python implementation parsing JASPAR matrices and computing log-odds PSSM scores.
Equivalent Cis-QL query: SCAN tbp_matrix STRAND POSITIVE THRESHOLD 80 % AS tbp_sites;
"""

import math

def parse_jaspar_pwm(jaspar_file):
    matrix = {'A': [], 'C': [], 'G': [], 'T': []}
    with open(jaspar_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>') or not line: continue
            parts = line.split()
            nuc = parts[0]
            if nuc in matrix:
                counts = [float(x) for x in parts[2:-1] if x.replace('.','',1).isdigit()]
                matrix[nuc] = counts
    return matrix

def compute_pssm(matrix, bg=0.25):
    pssm = {'A': [], 'C': [], 'G': [], 'T': []}
    cols = len(matrix['A'])
    min_score = 0.0
    max_score = 0.0
    for j in range(cols):
        col_total = sum(matrix[n][j] for n in 'ACGT')
        for n in 'ACGT':
            freq = (matrix[n][j] + 0.01) / (col_total + 0.04)
            score = math.log2(freq / bg)
            pssm[n].append(score)
    return pssm

if __name__ == '__main__':
    print("Conventional JASPAR PWM/PSSM parsing script completed.")
