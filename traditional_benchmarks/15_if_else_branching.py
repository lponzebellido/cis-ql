#!/usr/bin/env python3
"""
15_if_else_branching.py
Conventional Python script evaluating GC content and branching logic.
Equivalent Cis-QL Query: IF GC_CONTENT > 50 % THEN SCAN sp1 THRESHOLD 80 % AS gc_sites; ELSE SCAN tbp THRESHOLD 80 % AS at_sites; ENDIF;
"""

def if_else_branching(sequence):
    gc_count = sum(1 for b in sequence.upper() if b in 'GC')
    gc_pct = (gc_count / len(sequence)) * 100.0 if sequence else 0.0
    
    if gc_pct > 50.0:
        print(f"[Branch THEN] GC content is {gc_pct:.2f}% (> 50%). Scanning SP1 matrix...")
        selected_target = "SP1"
    else:
        print(f"[Branch ELSE] GC content is {gc_pct:.2f}% (<= 50%). Scanning TBP matrix...")
        selected_target = "TBP"
    return selected_target

if __name__ == '__main__':
    seq = "ATGCGCGCGCGCGCG"
    res = if_else_branching(seq)
