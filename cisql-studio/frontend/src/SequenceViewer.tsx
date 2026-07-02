import React from 'react';

interface GenomicRegion {
  chr: string;
  start: number;
  end: number;
  strand: string;
  type: string;
  name: string;
  sequence?: string;
}

interface SequenceViewerProps {
  results: Record<string, GenomicRegion[]>;
}

export const SequenceViewer: React.FC<SequenceViewerProps> = ({ results }) => {
  const keys = Object.keys(results);

  if (keys.length === 0) {
    return (
      <div className="welcome-screen" style={{ height: '100%', justifyContent: 'center', opacity: 0.5 }}>
        <p>No sequence results available.</p>
        <p>Run a query to populate this view.</p>
      </div>
    );
  }

  
  const MAX_DISPLAY = 500;
  
  return (
    <div className="sequence-viewer">
      {keys.map(track => (
        <div key={track} style={{ marginBottom: '32px' }}>
          <div style={{ color: 'var(--accent)', fontWeight: 600, borderBottom: '1px solid var(--border-color)', marginBottom: '16px', paddingBottom: '8px' }}>
            Result Set: {track} ({results[track].length} items)
          </div>
          
          {results[track].slice(0, MAX_DISPLAY).map((region, idx) => (
            <div key={idx} className="seq-region">
              <div className="seq-header">
                {region.chr} : {region.start}..{region.end} ({region.strand}) - {region.type}
              </div>
              <div style={{ wordBreak: 'break-all', fontFamily: "'JetBrains Mono', monospace", letterSpacing: '1px' }}>
                <span style={{ opacity: 0.5 }}>...</span>
                <span className="seq-highlight">
                  {region.sequence || 'N/A'}
                </span>
                <span style={{ opacity: 0.5 }}>...</span>
              </div>
            </div>
          ))}
          
          {results[track].length > MAX_DISPLAY && (
            <div style={{ color: 'var(--text-muted)', fontStyle: 'italic', marginTop: '16px' }}>
              Showing first {MAX_DISPLAY} results... ({results[track].length - MAX_DISPLAY} omitted for performance)
            </div>
          )}
        </div>
      ))}
    </div>
  );
};
