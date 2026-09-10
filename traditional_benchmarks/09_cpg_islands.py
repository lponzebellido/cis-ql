#!/usr/bin/env python3

def find_cpg_islands(sequence, window_size=200, step=100, min_gc=50.0, min_oe=0.6):
    islands = []
    n = len(sequence)
    for i in range(0, n - window_size, step):
        win = sequence[i:i+window_size].upper()
        g = win.count('G')
        c = win.count('C')
        cg = win.count('CG')
        gc_pct = ((g + c) / window_size) * 100.0
        exp_cg = (c * g) / window_size if window_size > 0 else 0
        obs_exp = (cg / exp_cg) if exp_cg > 0 else 0
        if gc_pct >= min_gc and obs_exp >= min_oe:
            islands.append((i, i + window_size, gc_pct, obs_exp))
    return islands

if __name__ == '__main__':
    seq = "CGCG" * 100
    islands = find_cpg_islands(seq)
    print(f"Discovered {len(islands)} CpG islands using conventional sliding window algorithm.")
