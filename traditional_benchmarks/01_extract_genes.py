#!/usr/bin/env python3
"""
01_extract_genes.py
Conventional Python implementation using imperative GFF3 file parsing and length filtering.
Equivalent Cis-QL query: EXTRACT GENE WHERE LENGTH > 1.5 KB;
"""

import sys

def extract_genes(gff_path, min_length_bp=1500):
    extracted_genes = []
    with open(gff_path, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            feature_type = parts[2]
            if feature_type == 'gene':
                start = int(parts[3])
                end = int(parts[4])
                length = end - start + 1
                if length > min_length_bp:
                    extracted_genes.append({
                        'chr': parts[0],
                        'start': start,
                        'end': end,
                        'length': length,
                        'strand': parts[6],
                        'attributes': parts[8]
                    })
    return extracted_genes

if __name__ == '__main__':
    genes = extract_genes("data_examples/genomic.gff", 1500)
    print(f"Extracted {len(genes)} structural genes (> 1.5 KB)")
