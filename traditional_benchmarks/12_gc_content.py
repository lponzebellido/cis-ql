#!/usr/bin/env python3
"""
12_gc_content.py
Conventional Python script calculating sliding-window GC content profiles.
Equivalent Cis-QL Query: ANALYZE GC_CONTENT WINDOW 1 KB AS gc_profile;
"""

def compute_gc_landscape(sequence, window_size=1000):
    landscape = []
    for i in range(0, len(sequence), window_size):
        chunk = sequence[i:i+window_size].upper()
        if not chunk: continue
        gc = sum(1 for b in chunk if b in 'GC')
        pct = (gc / len(chunk)) * 100.0
        landscape.append((i, pct))
    return landscape

if __name__ == '__main__':
    seq = "ATGC" * 1000
    landscape = compute_gc_landscape(seq, 1000)
    print(f"Computed {len(landscape)} GC windows across sequence.")
