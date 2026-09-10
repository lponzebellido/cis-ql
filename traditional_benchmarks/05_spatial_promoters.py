#!/usr/bin/env python3

import re

def spatial_search(gff_path, fasta_seq, motif="TTGACA", distance=200):
    cds_regions = []
    with open(gff_path, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip(): continue
            parts = line.strip().split('\t')
            if len(parts) >= 9 and parts[2] == 'CDS':
                cds_regions.append((int(parts[3]), int(parts[4]), parts[6]))
    
    spatial_hits = []
    for start, end, strand in cds_regions:
        if strand == '+':
            win_start = max(0, start - distance)
            win_end = start
        else:
            win_start = end
            win_end = min(len(fasta_seq), end + distance)
        
        subseq = fasta_seq[win_start:win_end]
        for m in re.finditer(motif, subseq):
            spatial_hits.append((win_start + m.start(), win_start + m.end()))
    return spatial_hits

if __name__ == '__main__':
    print("Conventional spatial promoter search completed.")
