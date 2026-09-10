#!/usr/bin/env python3

import re

def denovo_orf_discovery(fasta_seq, distance=300, min_orf_len=600):
    promoters = [m.end() for m in re.finditer(r"TATA[AT]A[AT]", fasta_seq)]
    candidate_orfs = []
    orf_regex = re.compile(r"ATG(?:...)*?(?:TAA|TAG|TGA)")
    
    for p_end in promoters:
        window_start = p_end
        window_end = min(len(fasta_seq), p_end + distance + 2000)
        subseq = fasta_seq[window_start:window_end]
        for match in orf_regex.finditer(subseq):
            orf_len = match.end() - match.start()
            if orf_len > min_orf_len:
                candidate_orfs.append((window_start + match.start(), window_start + match.end(), orf_len))
    return candidate_orfs

if __name__ == '__main__':
    seq = "TATAATAT" + "ATG" + "GCC" * 250 + "TAA"
    orfs = denovo_orf_discovery(seq)
    print(f"Discovered {len(orfs)} de novo ORFs using conventional Python regex pipeline.")
