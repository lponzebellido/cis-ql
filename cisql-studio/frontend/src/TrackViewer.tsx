import React, { useRef, useEffect, useState, useCallback } from 'react';

interface GenomicRegion {
  chr: string;
  start: number;
  end: number;
  strand: string;
  type: string;
  name: string;
  sequence?: string;
}

interface TrackViewerProps {
  results: Record<string, GenomicRegion[]>;
  onSelectRegion?: (region: GenomicRegion) => void;
}

const TYPE_COLORS: Record<string, string> = {
  gene: '#22863a',
  CDS: '#0366d6',
  exon: '#6f42c1',
  intron: '#959da5',
  UTR: '#e36209',
  promoter: '#d73a49',
  enhancer: '#b392f0',
  mRNA: '#28a745',
  region: '#f9c513',
  motif: '#ea4aaa',
};

function getTypeColor(type: string): string {
  const lower = type.toLowerCase();
  for (const [key, color] of Object.entries(TYPE_COLORS)) {
    if (lower.includes(key.toLowerCase())) return color;
  }
  const hash = Array.from(type).reduce((a, c) => a + c.charCodeAt(0), 0);
  const hue = hash % 360;
  return `hsl(${hue}, 60%, 55%)`;
}

function formatBp(bp: number): string {
  if (Math.abs(bp) >= 1e6) return (bp / 1e6).toFixed(2) + ' Mb';
  if (Math.abs(bp) >= 1e3) return (bp / 1e3).toFixed(1) + ' kb';
  return bp.toLocaleString() + ' bp';
}

export const TrackViewer: React.FC<TrackViewerProps> = ({ results, onSelectRegion }) => {
  const canvasRef = useRef<HTMLCanvasElement>(null);
  const containerRef = useRef<HTMLDivElement>(null);
  const [hoveredRegion, setHoveredRegion] = useState<GenomicRegion | null>(null);
  const [selectedRegion, setSelectedRegion] = useState<GenomicRegion | null>(null);
  const [mousePos, setMousePos] = useState({ x: 0, y: 0 });
  const [cursorBp, setCursorBp] = useState<number | null>(null);
  const [zoom, setZoom] = useState(1);
  const [panOffset, setPanOffset] = useState(0);
  const isDragging = useRef(false);
  const dragStartX = useRef(0);
  const panStartVal = useRef(0);

  const trackNames = React.useMemo(() => Object.keys(results), [results]);
  let globalMin = Infinity, globalMax = 0;
  trackNames.forEach(name => {
    results[name].forEach(r => {
      if (r.start < globalMin) globalMin = r.start;
      if (r.end > globalMax) globalMax = r.end;
    });
  });
  if (globalMin === Infinity) { globalMin = 0; globalMax = 1000; }
  const totalRange = globalMax - globalMin || 1000;
  const pad = totalRange * 0.05;

  const usedTypes = new Set<string>();
  trackNames.forEach(name => results[name].forEach(r => usedTypes.add(r.type)));

  const totalResults = trackNames.reduce((sum, n) => sum + results[n].length, 0);

  const RULER_H = 36;
  const LABEL_W = 90;
  const TRACK_H = 50;

  const draw = useCallback(() => {
    const canvas = canvasRef.current;
    const container = containerRef.current;
    if (!canvas || !container) return;
    const ctx = canvas.getContext('2d');
    if (!ctx) return;

    const rect = container.getBoundingClientRect();
    const dpr = window.devicePixelRatio || 1;
    canvas.width = rect.width * dpr;
    canvas.height = rect.height * dpr;
    ctx.scale(dpr, dpr);
    canvas.style.width = `${rect.width}px`;
    canvas.style.height = `${rect.height}px`;

    const w = rect.width;
    const h = rect.height;
    const dark = window.matchMedia('(prefers-color-scheme: dark)').matches;

    ctx.fillStyle = dark ? '#0d1117' : '#fff';
    ctx.fillRect(0, 0, w, h);

    if (trackNames.length === 0) return;

    const viewRange = totalRange / zoom;
    const viewMin = globalMin - pad + panOffset;
    const viewMax = viewMin + viewRange + pad * 2;
    const effRange = viewMax - viewMin;
    const toX = (bp: number) => LABEL_W + ((bp - viewMin) / effRange) * (w - LABEL_W);

    // Ruler background
    ctx.fillStyle = dark ? '#161b22' : '#f6f8fa';
    ctx.fillRect(0, 0, w, RULER_H);
    ctx.strokeStyle = dark ? '#30363d' : '#d0d7de';
    ctx.lineWidth = 1;
    ctx.beginPath();
    ctx.moveTo(LABEL_W, RULER_H - 0.5);
    ctx.lineTo(w, RULER_H - 0.5);
    ctx.stroke();

    // Ruler ticks
    const idealTicks = 8;
    const rawStep = effRange / idealTicks;
    const mag = Math.pow(10, Math.floor(Math.log10(rawStep)));
    const norm = rawStep / mag;
    let step: number;
    if (norm <= 1.5) step = 1 * mag;
    else if (norm <= 3) step = 2 * mag;
    else if (norm <= 7) step = 5 * mag;
    else step = 10 * mag;

    const firstTick = Math.ceil(viewMin / step) * step;
    ctx.textAlign = 'center';
    ctx.textBaseline = 'top';

    for (let bp = firstTick; bp <= viewMax; bp += step) {
      const x = toX(bp);
      if (x < LABEL_W || x > w) continue;

      ctx.strokeStyle = dark ? '#30363d' : '#d0d7de';
      ctx.beginPath();
      ctx.moveTo(x, RULER_H - 10);
      ctx.lineTo(x, RULER_H);
      ctx.stroke();

      ctx.fillStyle = dark ? '#8b949e' : '#57606a';
      ctx.font = '10px "JetBrains Mono", monospace';
      ctx.fillText(formatBp(bp), x, 4);
    }

    // Minor ticks
    const minorStep = step / 5;
    for (let bp = Math.ceil(viewMin / minorStep) * minorStep; bp <= viewMax; bp += minorStep) {
      const x = toX(bp);
      if (x < LABEL_W || x > w) continue;
      ctx.strokeStyle = dark ? '#21262d' : '#e8ebef';
      ctx.beginPath();
      ctx.moveTo(x, RULER_H - 4);
      ctx.lineTo(x, RULER_H);
      ctx.stroke();
    }

    // Label column background
    ctx.fillStyle = dark ? '#161b22' : '#f6f8fa';
    ctx.fillRect(0, RULER_H, LABEL_W, h - RULER_H);
    ctx.strokeStyle = dark ? '#30363d' : '#d0d7de';
    ctx.beginPath();
    ctx.moveTo(LABEL_W - 0.5, RULER_H);
    ctx.lineTo(LABEL_W - 0.5, h);
    ctx.stroke();

    // Tracks
    trackNames.forEach((trackName, idx) => {
      const yTop = RULER_H + idx * TRACK_H;
      const yMid = yTop + TRACK_H / 2;

      // Track separator
      if (idx > 0) {
        ctx.strokeStyle = dark ? '#21262d' : '#eaeef2';
        ctx.beginPath();
        ctx.moveTo(0, yTop);
        ctx.lineTo(w, yTop);
        ctx.stroke();
      }

      // Alternating track bg
      if (idx % 2 === 1) {
        ctx.fillStyle = dark ? 'rgba(22,27,34,0.5)' : 'rgba(246,248,250,0.5)';
        ctx.fillRect(LABEL_W, yTop, w - LABEL_W, TRACK_H);
      }

      // Track label
      ctx.fillStyle = dark ? '#e6edf3' : '#1f2328';
      ctx.font = '11px Inter, sans-serif';
      ctx.textAlign = 'left';
      ctx.textBaseline = 'middle';
      ctx.fillText(trackName.length > 10 ? trackName.substring(0, 10) + '..' : trackName, 8, yMid);

      // Center line
      ctx.strokeStyle = dark ? '#21262d' : '#e1e4e8';
      ctx.setLineDash([3, 3]);
      ctx.beginPath();
      ctx.moveTo(LABEL_W, yMid);
      ctx.lineTo(w, yMid);
      ctx.stroke();
      ctx.setLineDash([]);

      // Features
      const featureH = TRACK_H * 0.32;
      results[trackName].forEach(r => {
        const x1 = toX(r.start);
        const x2 = toX(r.end);
        if (x2 < LABEL_W || x1 > w) return;
        const clampX1 = Math.max(LABEL_W, x1);
        const fW = Math.max(2, x2 - clampX1);
        const color = getTypeColor(r.type);
        const yOff = r.strand === '+' ? yMid - featureH - 2 : yMid + 2;

        const isSelected = selectedRegion && selectedRegion.name === r.name && selectedRegion.start === r.start;
        const isHovered = hoveredRegion && hoveredRegion.name === r.name && hoveredRegion.start === r.start;

        ctx.fillStyle = color;
        ctx.globalAlpha = isSelected ? 1.0 : isHovered ? 0.9 : 0.75;

        if (fW > 12) {
          const arrowW = Math.min(8, fW * 0.12);
          ctx.beginPath();
          if (r.strand === '+') {
            ctx.moveTo(clampX1, yOff);
            ctx.lineTo(clampX1 + fW - arrowW, yOff);
            ctx.lineTo(clampX1 + fW, yOff + featureH / 2);
            ctx.lineTo(clampX1 + fW - arrowW, yOff + featureH);
            ctx.lineTo(clampX1, yOff + featureH);
          } else {
            ctx.moveTo(clampX1 + arrowW, yOff);
            ctx.lineTo(clampX1 + fW, yOff);
            ctx.lineTo(clampX1 + fW, yOff + featureH);
            ctx.lineTo(clampX1 + arrowW, yOff + featureH);
            ctx.lineTo(clampX1, yOff + featureH / 2);
          }
          ctx.closePath();
          ctx.fill();

          if (isSelected) {
            ctx.strokeStyle = '#fff';
            ctx.lineWidth = 2;
            ctx.stroke();
            ctx.lineWidth = 1;
          }
        } else {
          ctx.fillRect(clampX1, yOff, fW, featureH);
          if (isSelected) {
            ctx.strokeStyle = '#fff';
            ctx.lineWidth = 2;
            ctx.strokeRect(clampX1, yOff, fW, featureH);
            ctx.lineWidth = 1;
          }
        }

        ctx.globalAlpha = 1.0;

        // Label inside feature if wide enough
        if (fW > 50) {
          ctx.fillStyle = '#fff';
          ctx.font = '9px Inter, sans-serif';
          ctx.textAlign = 'center';
          ctx.textBaseline = 'middle';
          const label = r.name.length > Math.floor(fW / 7) ? r.name.substring(0, Math.floor(fW / 7)) + '..' : r.name;
          ctx.fillText(label, clampX1 + fW / 2, yOff + featureH / 2);
        }
      });
    });

    // Cursor crosshair
    if (cursorBp !== null) {
      const cx = toX(cursorBp);
      if (cx > LABEL_W && cx < w) {
        ctx.strokeStyle = dark ? 'rgba(88,166,255,0.4)' : 'rgba(9,105,218,0.3)';
        ctx.setLineDash([4, 4]);
        ctx.beginPath();
        ctx.moveTo(cx, RULER_H);
        ctx.lineTo(cx, h);
        ctx.stroke();
        ctx.setLineDash([]);
      }
    }
  }, [results, zoom, panOffset, trackNames, globalMin, globalMax, totalRange, pad, hoveredRegion, selectedRegion, cursorBp]);

  useEffect(() => { draw(); }, [draw]);

  useEffect(() => {
    const obs = new ResizeObserver(() => draw());
    if (containerRef.current) obs.observe(containerRef.current);
    return () => obs.disconnect();
  }, [draw]);

  const getViewParams = () => {
    const rect = containerRef.current?.getBoundingClientRect();
    if (!rect) return null;
    const viewRange = totalRange / zoom;
    const viewMin = globalMin - pad + panOffset;
    const viewMax = viewMin + viewRange + pad * 2;
    const effRange = viewMax - viewMin;
    return { viewMin, effRange, w: rect.width, h: rect.height };
  };

  const handleWheel = (e: React.WheelEvent) => {
    e.preventDefault();
    const factor = e.deltaY > 0 ? 0.85 : 1.18;
    setZoom(prev => Math.max(0.5, Math.min(200, prev * factor)));
  };

  const handleMouseDown = (e: React.MouseEvent) => {
    isDragging.current = true;
    dragStartX.current = e.clientX;
    panStartVal.current = panOffset;
    (e.target as HTMLElement).style.cursor = 'grabbing';
  };

  const handleMouseMove = (e: React.MouseEvent<HTMLCanvasElement>) => {
    const p = getViewParams();
    if (!p) return;
    const rect = containerRef.current!.getBoundingClientRect();
    const x = e.clientX - rect.left;
    const y = e.clientY - rect.top;

    const bp = p.viewMin + ((x - LABEL_W) / (p.w - LABEL_W)) * p.effRange;
    setCursorBp(bp);

    if (isDragging.current) {
      const dx = e.clientX - dragStartX.current;
      const bpPerPx = p.effRange / (p.w - LABEL_W);
      setPanOffset(panStartVal.current - dx * bpPerPx);
      return;
    }

    const trackIdx = Math.floor((y - RULER_H) / TRACK_H);
    let found: GenomicRegion | null = null;
    if (trackIdx >= 0 && trackIdx < trackNames.length && x > LABEL_W) {
      const hitPad = p.effRange * 0.002;
      for (const r of results[trackNames[trackIdx]]) {
        if (bp >= r.start - hitPad && bp <= r.end + hitPad) {
          found = r;
          break;
        }
      }
    }
    setHoveredRegion(found);
    setMousePos({ x: e.clientX, y: e.clientY });
  };

  const handleMouseUp = (e: React.MouseEvent) => {
    const wasDragging = isDragging.current && Math.abs(e.clientX - dragStartX.current) > 3;
    isDragging.current = false;
    (e.target as HTMLElement).style.cursor = '';

    if (!wasDragging && hoveredRegion) {
      setSelectedRegion(hoveredRegion);
      if (onSelectRegion) onSelectRegion(hoveredRegion);
    }
  };

  const handleMouseLeave = () => {
    setHoveredRegion(null);
    setCursorBp(null);
    isDragging.current = false;
  };

  const zoomToRegion = (region: GenomicRegion) => {
    const regionRange = region.end - region.start;
    const buffer = regionRange * 2;
    const center = (region.start + region.end) / 2;
    const newZoom = totalRange / (regionRange + buffer);
    setZoom(Math.min(200, Math.max(0.5, newZoom)));
    setPanOffset(center - globalMin - totalRange / (2 * newZoom));
  };

  if (trackNames.length === 0) {
    return (
      <div className="track-viewer-container">
        <div className="track-empty">
          <div>
            <div style={{ fontSize: '1.1rem', marginBottom: 8 }}>No regions to display</div>
            <div style={{ fontSize: '0.8rem', opacity: 0.6 }}>Run a Cis-QL query to visualize genomic results here.</div>
          </div>
        </div>
      </div>
    );
  }

  return (
    <div className="track-viewer-container">
      <div className="track-toolbar">
        <button onClick={() => setZoom(prev => Math.min(200, prev * 1.5))}>+ Zoom In</button>
        <button onClick={() => setZoom(prev => Math.max(0.5, prev / 1.5))}>- Zoom Out</button>
        <button onClick={() => { setZoom(1); setPanOffset(0); setSelectedRegion(null); }}>Reset</button>
        {selectedRegion && (
          <button onClick={() => zoomToRegion(selectedRegion)}>Focus Selected</button>
        )}
        <span className="zoom-info">
          {zoom.toFixed(1)}x | {totalResults} features | {trackNames.length} tracks
          {cursorBp !== null && ` | ${formatBp(Math.round(cursorBp))}`}
        </span>
      </div>

      <div className="track-legend">
        {Array.from(usedTypes).map(type => (
          <div key={type} className="legend-item">
            <div className="legend-color" style={{ background: getTypeColor(type) }}></div>
            {type}
          </div>
        ))}
        <div className="legend-item" style={{ marginLeft: 'auto', opacity: 0.6 }}>
          Top = Forward (+) | Bottom = Reverse (-)
        </div>
      </div>

      <div style={{ display: 'flex', flex: 1, minHeight: 0 }}>
        <div className="track-canvas-wrap" ref={containerRef}>
          <canvas
            ref={canvasRef}
            onWheel={handleWheel}
            onMouseDown={handleMouseDown}
            onMouseMove={handleMouseMove}
            onMouseUp={handleMouseUp}
            onMouseLeave={handleMouseLeave}
            style={{ cursor: hoveredRegion ? 'pointer' : 'grab' }}
          />
          {hoveredRegion && !isDragging.current && (
            <div className="track-tooltip" style={{ left: mousePos.x + 14, top: mousePos.y - 10 }}>
              <strong>{hoveredRegion.name}</strong><br />
              <span className="tt-label">Type:</span> {hoveredRegion.type}<br />
              <span className="tt-label">Pos:</span> {hoveredRegion.chr}:{hoveredRegion.start.toLocaleString()}-{hoveredRegion.end.toLocaleString()}<br />
              <span className="tt-label">Strand:</span> {hoveredRegion.strand === '+' ? 'Forward (+)' : 'Reverse (-)'}<br />
              <span className="tt-label">Length:</span> {(hoveredRegion.end - hoveredRegion.start).toLocaleString()} bp
              {hoveredRegion.sequence && (
                <div style={{ marginTop: 6, fontFamily: "'JetBrains Mono', monospace", fontSize: '0.65rem', color: '#3fb950', wordBreak: 'break-all', maxWidth: 250 }}>
                  {hoveredRegion.sequence.substring(0, 40)}{hoveredRegion.sequence.length > 40 ? '...' : ''}
                </div>
              )}
              <div style={{ marginTop: 4, fontSize: '0.6rem', opacity: 0.5 }}>Click to select</div>
            </div>
          )}
        </div>

        {selectedRegion && (
          <div className="detail-panel">
            <div className="detail-header">
              <span>Selected Region</span>
              <button onClick={() => setSelectedRegion(null)} style={{ background: 'none', border: 'none', color: 'var(--text-muted)', cursor: 'pointer' }}>x</button>
            </div>
            <div className="detail-body">
              <div className="detail-field">
                <span className="detail-label">Name</span>
                <span className="detail-value">{selectedRegion.name}</span>
              </div>
              <div className="detail-field">
                <span className="detail-label">Type</span>
                <span className="detail-value" style={{ display: 'flex', alignItems: 'center', gap: 6 }}>
                  <span style={{ width: 8, height: 8, borderRadius: 2, background: getTypeColor(selectedRegion.type), display: 'inline-block' }}></span>
                  {selectedRegion.type}
                </span>
              </div>
              <div className="detail-field">
                <span className="detail-label">Chromosome</span>
                <span className="detail-value">{selectedRegion.chr}</span>
              </div>
              <div className="detail-field">
                <span className="detail-label">Position</span>
                <span className="detail-value">{selectedRegion.start.toLocaleString()} - {selectedRegion.end.toLocaleString()}</span>
              </div>
              <div className="detail-field">
                <span className="detail-label">Length</span>
                <span className="detail-value">{formatBp(selectedRegion.end - selectedRegion.start)}</span>
              </div>
              <div className="detail-field">
                <span className="detail-label">Strand</span>
                <span className="detail-value">{selectedRegion.strand === '+' ? 'Forward (+)' : 'Reverse (-)'}</span>
              </div>
              {selectedRegion.sequence && (
                <div className="detail-seq">
                  <div className="detail-label" style={{ marginBottom: 4 }}>Sequence ({selectedRegion.sequence.length} bp)</div>
                  <div className="detail-seq-text" onClick={() => {
                    navigator.clipboard.writeText(selectedRegion.sequence || '');
                  }} title="Click to copy">
                    {selectedRegion.sequence}
                  </div>
                  <div style={{ fontSize: '0.6rem', color: 'var(--text-muted)', marginTop: 4 }}>Click sequence to copy to clipboard</div>
                </div>
              )}
              <button className="detail-action" onClick={() => zoomToRegion(selectedRegion)}>
                Zoom to Region
              </button>
            </div>
          </div>
        )}
      </div>
    </div>
  );
};
