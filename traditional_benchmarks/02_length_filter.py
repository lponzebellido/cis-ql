#!/usr/bin/env python3

def filter_genes_by_range(gff_path, min_len=500, max_len=3000):
    results = []
    with open(gff_path, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 9 and parts[2] == 'gene':
                start = int(parts[3])
                end = int(parts[4])
                length = end - start + 1
                if min_len <= length <= max_len:
                    results.append((parts[0], start, end, length))
    return results

if __name__ == '__main__':
    res = filter_genes_by_range("data_examples/genomic.gff", 500, 3000)
    print(f"Found {len(res)} genes within range [500 BP, 3 KB]")
