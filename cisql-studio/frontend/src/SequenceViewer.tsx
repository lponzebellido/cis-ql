import React, { useState, useMemo } from 'react';

interface MotifEvidence {
  matrixAlias: string;
  matrixId: string;
  matrixName: string;
  matrixSource: string;
  rawScore: number;
  scorePercent: number;
  sourceRegion?: {
    name: string;
    type: string;
    start: number;
    end: number;
    relativeStart: number;
  };
}

interface GenomicRegion {
  chr: string;
  start: number;
  end: number;
  strand: string;
  type: string;
  name: string;
  sequence?: string;
  motifEvidence?: MotifEvidence;
}

interface SequenceViewerProps {
  results: Record<string, GenomicRegion[]>;
  highlightedRegion?: GenomicRegion | null;
}

function ntClass(base: string): string {
  switch (base.toUpperCase()) {
    case 'A': return 'nt-A';
    case 'T': return 'nt-T';
    case 'G': return 'nt-G';
    case 'C': return 'nt-C';
    default: return 'nt-other';
  }
}

function gcContent(seq: string): string {
  if (!seq) return 'N/A';
  let gc = 0;
  for (const c of seq) {
    if (c === 'G' || c === 'C' || c === 'g' || c === 'c') gc++;
  }
  return ((gc / seq.length) * 100).toFixed(1) + '%';
}

function ntCounts(seq: string): Record<string, number> {
  const counts: Record<string, number> = { A: 0, T: 0, G: 0, C: 0, Other: 0 };
  for (const c of seq) {
    const upper = c.toUpperCase();
    if (upper in counts) counts[upper]++;
    else counts['Other']++;
  }
  return counts;
}

function renderSequenceBlock(seq: string, startPos: number) {
  const BASES_PER_LINE = 60;
  const BLOCK_SIZE = 10;
  const lines: React.ReactNode[] = [];

  for (let i = 0; i < seq.length; i += BASES_PER_LINE) {
    const lineSeq = seq.substring(i, i + BASES_PER_LINE);
    const pos = startPos + i;
    const blocks: React.ReactNode[] = [];

    for (let b = 0; b < lineSeq.length; b += BLOCK_SIZE) {
      const blockStr = lineSeq.substring(b, b + BLOCK_SIZE);
      blocks.push(
        <span key={b} className="seq-block">
          {blockStr.split('').map((base, j) => (
            <span key={j} className={ntClass(base)}>{base}</span>
          ))}
        </span>
      );
    }

    lines.push(
      <div key={i} className="seq-line">
        <span className="seq-pos">{pos.toLocaleString()}</span>
        <span className="seq-bases">{blocks}</span>
        <span className="seq-pos-end">{Math.min(pos + BASES_PER_LINE - 1, startPos + seq.length - 1).toLocaleString()}</span>
      </div>
    );
  }
  return lines;
}

export const SequenceViewer: React.FC<SequenceViewerProps> = ({ results, highlightedRegion }) => {
  const [searchQuery, setSearchQuery] = useState('');
  const [expandedIdx, setExpandedIdx] = useState<string | null>(null);
  const [copiedIdx, setCopiedIdx] = useState<string | null>(null);

  const keys = Object.keys(results);

  const stats = useMemo(() => {
    const s: Record<string, { total: number; withSeq: number; avgLen: number; types: Set<string> }> = {};
    keys.forEach(track => {
      const regions = results[track];
      let seqCount = 0;
      let totalLen = 0;
      const types = new Set<string>();
      regions.forEach(r => {
        if (r.sequence) seqCount++;
        totalLen += r.end - r.start;
        types.add(r.type);
      });
      s[track] = {
        total: regions.length,
        withSeq: seqCount,
        avgLen: regions.length > 0 ? Math.round(totalLen / regions.length) : 0,
        types,
      };
    });
    return s;
  }, [results, keys]);

  const filteredResults = useMemo(() => {
    if (!searchQuery.trim()) return results;
    const q = searchQuery.toLowerCase();
    const filtered: Record<string, GenomicRegion[]> = {};
    keys.forEach(track => {
      filtered[track] = results[track].filter(r =>
        r.name.toLowerCase().includes(q) ||
        r.type.toLowerCase().includes(q) ||
        r.chr.toLowerCase().includes(q) ||
        (r.sequence && r.sequence.toUpperCase().includes(q.toUpperCase()))
      );
    });
    return filtered;
  }, [results, keys, searchQuery]);

  const handleCopy = async (seq: string, idx: string) => {
    try {
      await navigator.clipboard.writeText(seq);
      setCopiedIdx(idx);
      setTimeout(() => setCopiedIdx(null), 1500);
    } catch {
    }
  };

  const handleExportFasta = (track: string) => {
    const regions = results[track].filter(r => r.sequence);
    if (regions.length === 0) return;
    let fasta = '';
    regions.forEach(r => {
      fasta += `>${r.name} ${r.chr}:${r.start}-${r.end} (${r.strand}) ${r.type}\n`;
      const seq = r.sequence || '';
      for (let i = 0; i < seq.length; i += 80) {
        fasta += seq.substring(i, i + 80) + '\n';
      }
    });
    const blob = new Blob([fasta], { type: 'text/plain' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = `${track}_results.fasta`;
    a.click();
    URL.revokeObjectURL(url);
  };

  if (keys.length === 0) {
    return (
      <div className="seq-empty">
        <div style={{ fontSize: '1.1rem', marginBottom: 8 }}>No sequence results available</div>
        <div style={{ fontSize: '0.8rem', opacity: 0.6 }}>Run a Cis-QL query to populate this view.</div>
      </div>
    );
  }

  const MAX_DISPLAY = 100;

  return (
    <div className="sequence-viewer">
      <div className="seq-toolbar">
        <input
          className="seq-search"
          type="text"
          placeholder="Filter by name, type, chr, or sequence..."
          value={searchQuery}
          onChange={e => setSearchQuery(e.target.value)}
        />
      </div>

      {keys.map(track => {
        const trackStats = stats[track];
        const regions = filteredResults[track] || [];
        return (
          <div key={track} className="seq-track-group">
            <div className="seq-track-header">
              <div className="seq-track-title">
                <span className="seq-track-name">{track}</span>
                <span className="seq-track-count">{regions.length} results</span>
              </div>
              <div className="seq-track-stats">
                <span>Avg length: {trackStats.avgLen.toLocaleString()} bp</span>
                <span>With sequence: {trackStats.withSeq}/{trackStats.total}</span>
                <span>Types: {Array.from(trackStats.types).join(', ')}</span>
              </div>
              <div className="seq-track-actions">
                {trackStats.withSeq > 0 && (
                  <button className="seq-action-btn" onClick={() => handleExportFasta(track)}>
                    Export FASTA
                  </button>
                )}
              </div>
            </div>

            <div className="seq-results-list">
              {regions.slice(0, MAX_DISPLAY).map((region, idx) => {
                const uid = `${track}-${idx}`;
                const isExpanded = expandedIdx === uid;
                const isHighlighted = highlightedRegion &&
                  highlightedRegion.name === region.name &&
                  highlightedRegion.start === region.start;

                return (
                  <div key={uid}
                    className={`seq-result-card ${isHighlighted ? 'highlighted' : ''}`}
                    id={isHighlighted ? 'highlighted-region' : undefined}
                  >
                    <div className="seq-result-header" onClick={() => setExpandedIdx(isExpanded ? null : uid)}>
                      <span className="seq-expand-icon">{isExpanded ? '▾' : '▸'}</span>
                      <span className="seq-result-name">{region.name}</span>
                      <span className="seq-result-location">
                        {region.chr}:{region.start.toLocaleString()}-{region.end.toLocaleString()}
                      </span>
                      <span className={`strand-badge ${region.strand === '+' ? 'forward' : 'reverse'}`}>
                        {region.strand === '+' ? 'Fwd' : 'Rev'}
                      </span>
                      <span className="seq-result-type">{region.type}</span>
                      <span className="seq-result-len">{(region.end - region.start).toLocaleString()} bp</span>
                    </div>

                    {isExpanded && (
                      <div className="seq-result-body">
                        {region.sequence ? (
                          <>
                            <div className="seq-body-stats">
                              <span>Length: {region.sequence.length} bp</span>
                              <span>GC Content: {gcContent(region.sequence)}</span>
                              {region.motifEvidence && (
                                <>
                                  <span>Matrix: {region.motifEvidence.matrixId || region.motifEvidence.matrixAlias}</span>
                                  <span>PWM score: {region.motifEvidence.rawScore.toFixed(3)} ({region.motifEvidence.scorePercent.toFixed(1)}%)</span>
                                  {region.motifEvidence.sourceRegion && (
                                    <span>Source: {region.motifEvidence.sourceRegion.name} +{region.motifEvidence.sourceRegion.relativeStart} bp</span>
                                  )}
                                </>
                              )}
                              {(() => {
                                const c = ntCounts(region.sequence!);
                                return (
                                  <>
                                    <span className="nt-A">A: {c.A}</span>
                                    <span className="nt-T">T: {c.T}</span>
                                    <span className="nt-G">G: {c.G}</span>
                                    <span className="nt-C">C: {c.C}</span>
                                  </>
                                );
                              })()}
                              <button
                                className="seq-copy-btn"
                                onClick={() => handleCopy(region.sequence!, uid)}
                              >
                                {copiedIdx === uid ? 'Copied' : 'Copy Sequence'}
                              </button>
                            </div>
                            <div className="seq-body-sequence">
                              {renderSequenceBlock(region.sequence, region.start)}
                            </div>
                          </>
                        ) : (
                          <div className="seq-no-data">No sequence data for this region.</div>
                        )}
                      </div>
                    )}
                  </div>
                );
              })}
            </div>

            {regions.length > MAX_DISPLAY && (
              <div className="seq-truncated">
                Showing {MAX_DISPLAY} of {regions.length} results.
              </div>
            )}
          </div>
        );
      })}
    </div>
  );
};
