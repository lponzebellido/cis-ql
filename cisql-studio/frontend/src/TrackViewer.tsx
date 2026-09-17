import React, { useRef, useEffect, useState, useCallback } from 'react';

interface MotifEvidence {
  matrixAlias: string;
  matrixId: string;
  matrixName: string;
  matrixSource: string;
  rawScore: number;
  scorePercent: number;
  motifPseudocount: number;
  statistics?: {
    pValue: number;
    qValue: number;
    testedPositions: number;
    pValueMethod: string;
    multipleTestingMethod: string;
    scaledScore: number;
    scoreRange: number;
    scoreScale: number;
    scoreOffset: number;
  };
  background?: {
    mode: string;
    source: string;
    A: number;
    C: number;
    G: number;
    T: number;
    estimationPseudocount: number;
    observedBases: number;
    strandPolicy: string;
  };
  sourceRegion?: {
    name: string;
    type: string;
    start: number;
    end: number;
    relativeStart: number;
  };
}

interface SpatialRelation {
  relation: string;
  referenceSet: string;
  reference: {
    chr: string;
    start: number;
    end: number;
    strand: string;
    type: string;
    name: string;
  };
  distance: number;
  maximumDistance: number;
  overlaps: boolean;
}

interface CountEvidence {
  relation: string;
  countedSet: string;
  containerSet: string;
  count: number;
}

interface ConsensusEvidence {
  anchorSet: string;
  minimumSupport: number;
  observedSupport: number;
  minimumReciprocalOverlapPercent?: number;
  inputSets: string[];
}

interface TrackEvidence {
  trackAlias: string;
  source: string;
  format: string;
  evidenceClass: string;
  assay?: string;
  sample?: string;
  condition?: string;
  replicate?: string;
  control?: string;
  score?: number;
  signalValue?: number;
  minusLog10PValue?: number;
  minusLog10QValue?: number;
  peakOffset?: number;
  peakPosition?: number;
}

interface OverlapEvidence {
  referenceSet: string;
  reference: {
    chr: string;
    start: number;
    end: number;
    strand: string;
    type: string;
    name: string;
  };
  trackEvidence?: TrackEvidence;
  consensusEvidence?: ConsensusEvidence;
  supportingEvidence?: OverlapEvidence[];
}

interface ModuleMemberEvidence {
  sourceSet: string;
  chr: string;
  start: number;
  end: number;
  strand: string;
  type: string;
  name: string;
  motifEvidence?: MotifEvidence;
}

interface ModuleEvidence {
  spacing: { minimum: number; maximum: number; observed: number };
  order: { policy: string; observed: string };
  orientation: { policy: string; observed: string };
  members: ModuleMemberEvidence[];
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
  spatialRelation?: SpatialRelation;
  countEvidence?: CountEvidence;
  consensusEvidence?: ConsensusEvidence;
  moduleEvidence?: ModuleEvidence;
  trackEvidence?: TrackEvidence;
  overlapEvidence?: OverlapEvidence[];
}

interface TrackViewerProps {
  results: Record<string, GenomicRegion[]>;
  gcProfiles?: Record<string, any[]>;
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
  cpg_island: '#00fa9a',
  gc_content: '#a2ee31',
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

function formatStrand(strand: string): string {
  if (strand === '+') return 'Forward (+)';
  if (strand === '-') return 'Reverse (-)';
  return 'Unstranded (.)';
}

export const TrackViewer: React.FC<TrackViewerProps> = ({ results, gcProfiles = {}, onSelectRegion }) => {
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
  const gcProfileNames = React.useMemo(() => Object.keys(gcProfiles), [gcProfiles]);
  let globalMin = Infinity, globalMax = 0;
  trackNames.forEach(name => {
    results[name].forEach(r => {
      if (r.start < globalMin) globalMin = r.start;
      if (r.end > globalMax) globalMax = r.end;
    });
  });
  gcProfileNames.forEach(name => {
    if (gcProfiles[name].length > 0) {
      if (gcProfiles[name][0].pos < globalMin) globalMin = gcProfiles[name][0].pos;
      if (gcProfiles[name][gcProfiles[name].length - 1].pos > globalMax) globalMax = gcProfiles[name][gcProfiles[name].length - 1].pos;
    }
  });

  if (globalMin === Infinity) { globalMin = 0; globalMax = 1000; }
  const totalRange = globalMax - globalMin || 1000;
  const pad = totalRange * 0.05;

  const usedTypes = new Set<string>();
  trackNames.forEach(name => results[name].forEach(r => usedTypes.add(r.type)));
  if (gcProfileNames.length > 0) usedTypes.add("gc_content");

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
    const requiredH = RULER_H + (gcProfileNames.length * 60) + (trackNames.length * 50) + 30;
    const contentH = Math.max(rect.height, requiredH);

    canvas.width = rect.width * dpr;
    canvas.height = contentH * dpr;
    ctx.scale(dpr, dpr);
    canvas.style.width = `${rect.width}px`;
    canvas.style.height = `${contentH}px`;

    const w = rect.width;
    const h = contentH;
    const dark = window.matchMedia('(prefers-color-scheme: dark)').matches;

    ctx.fillStyle = dark ? '#0d1117' : '#fff';
    ctx.fillRect(0, 0, w, h);

    if (trackNames.length === 0 && gcProfileNames.length === 0) return;

    const viewRange = totalRange / zoom;
    const padVal = viewRange * 0.02;
    const viewMin = globalMin - padVal + panOffset;
    const viewMax = viewMin + viewRange + padVal * 2;
    const effRange = viewMax - viewMin;
    const toX = (bp: number) => LABEL_W + ((bp - viewMin) / effRange) * (w - LABEL_W);

    ctx.fillStyle = dark ? '#161b22' : '#f6f8fa';
    ctx.fillRect(0, 0, w, RULER_H);
    ctx.strokeStyle = dark ? '#30363d' : '#d0d7de';
    ctx.lineWidth = 1;
    ctx.beginPath();
    ctx.moveTo(LABEL_W, RULER_H - 0.5);
    ctx.lineTo(w, RULER_H - 0.5);
    ctx.stroke();

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

    ctx.fillStyle = dark ? '#161b22' : '#f6f8fa';
    ctx.fillRect(0, RULER_H, LABEL_W, h - RULER_H);
    ctx.strokeStyle = dark ? '#30363d' : '#d0d7de';
    ctx.beginPath();
    ctx.moveTo(LABEL_W - 0.5, RULER_H);
    ctx.lineTo(LABEL_W - 0.5, h);
    ctx.stroke();

    let currentY = RULER_H;
    const GC_TRACK_H = 60;

    gcProfileNames.forEach((profileName, idx) => {
      const yTop = currentY;
      if (idx > 0 || trackNames.length > 0) {
        ctx.strokeStyle = dark ? '#21262d' : '#eaeef2';
        ctx.beginPath();
        ctx.moveTo(0, yTop);
        ctx.lineTo(w, yTop);
        ctx.stroke();
      }

      ctx.fillStyle = dark ? 'rgba(22,27,34,0.3)' : 'rgba(246,248,250,0.3)';
      ctx.fillRect(LABEL_W, yTop, w - LABEL_W, GC_TRACK_H);

      ctx.fillStyle = dark ? '#e6edf3' : '#1f2328';
      ctx.font = 'bold 11px Inter, sans-serif';
      ctx.textAlign = 'left';
      ctx.textBaseline = 'middle';
      ctx.fillText(profileName.length > 10 ? profileName.substring(0, 10) + '..' : profileName, 8, yTop + 15);
      ctx.fillStyle = dark ? '#8b949e' : '#57606a';
      ctx.font = '9px Inter, sans-serif';
      ctx.fillText('GC %', 8, yTop + 30);

      const profile = gcProfiles[profileName];
      if (profile.length > 0) {
        ctx.strokeStyle = getTypeColor('gc_content');
        ctx.fillStyle = dark ? 'rgba(162, 238, 49, 0.2)' : 'rgba(162, 238, 49, 0.3)';
        ctx.lineWidth = 1.5;

        ctx.beginPath();
        let started = false;

        for (let i = 0; i < profile.length; i++) {
          const w = profile[i];
          const x = toX(w.pos);
          if (x < LABEL_W) continue;
          if (!started) {
            ctx.moveTo(x, yTop + GC_TRACK_H - (w.gc / 100) * GC_TRACK_H);
            started = true;
          } else {
            ctx.lineTo(x, yTop + GC_TRACK_H - (w.gc / 100) * GC_TRACK_H);
          }
          if (x > w + 50) break;
        }
        ctx.stroke();

        if (started) {
          ctx.lineTo(toX(profile[profile.length - 1].pos), yTop + GC_TRACK_H);
          ctx.lineTo(toX(profile[0].pos), yTop + GC_TRACK_H);
          ctx.closePath();
          ctx.fill();
        }
      }
      currentY += GC_TRACK_H;
    });

    trackNames.forEach((trackName, idx) => {
      const yTop = currentY;
      const yMid = yTop + TRACK_H / 2;

      if (idx > 0) {
        ctx.strokeStyle = dark ? '#21262d' : '#eaeef2';
        ctx.beginPath();
        ctx.moveTo(0, yTop);
        ctx.lineTo(w, yTop);
        ctx.stroke();
      }

      if (idx % 2 === 1) {
        ctx.fillStyle = dark ? 'rgba(22,27,34,0.5)' : 'rgba(246,248,250,0.5)';
        ctx.fillRect(LABEL_W, yTop, w - LABEL_W, TRACK_H);
      }

      const trackFeats = results[trackName];
      let sumLen = 0;
      trackFeats.forEach(f => sumLen += (f.end - f.start));
      const avgLen = trackFeats.length ? Math.round(sumLen / trackFeats.length) : 0;
      ctx.fillStyle = dark ? '#e6edf3' : '#1f2328';
      ctx.font = 'bold 11px Inter, sans-serif';
      ctx.textAlign = 'left';
      ctx.textBaseline = 'middle';
      ctx.fillText(trackName.length > 10 ? trackName.substring(0, 10) + '..' : trackName, 8, yTop + 12);
      ctx.fillStyle = dark ? '#8b949e' : '#57606a';
      ctx.font = '9px Inter, sans-serif';
      ctx.fillText(`${trackFeats.length} features`, 8, yTop + 26);
      ctx.fillText(`~${formatBp(avgLen)} avg`, 8, yTop + 38);

      ctx.strokeStyle = dark ? '#21262d' : '#e1e4e8';
      ctx.setLineDash([3, 3]);
      ctx.beginPath();
      ctx.moveTo(LABEL_W, yMid);
      ctx.lineTo(w, yMid);
      ctx.stroke();
      ctx.setLineDash([]);

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
        } else if (fW < 1) {
          ctx.globalAlpha = isSelected ? 1.0 : isHovered ? 0.9 : 0.4;
          ctx.fillRect(clampX1, yOff, 1, featureH);
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

        const pxPerBp = fW / (r.end - r.start);
        if (pxPerBp > 8 && r.sequence) {
          ctx.font = 'bold 10px "JetBrains Mono", monospace';
          ctx.textAlign = 'center';
          ctx.textBaseline = 'middle';

          for (let i = 0; i < r.sequence.length; i++) {
            const base = r.sequence[i].toUpperCase();
            const baseX = toX(r.start + i) + pxPerBp / 2;
            if (baseX < LABEL_W || baseX > w) continue;

            ctx.fillStyle = base === 'A' ? '#3fb950' :
              base === 'T' ? '#f85149' :
                base === 'G' ? '#d29922' :
                  base === 'C' ? '#58a6ff' : '#8b949e';
            ctx.fillText(base, baseX, yOff + featureH / 2);
          }
        } else if (fW > 50) {
          ctx.fillStyle = '#fff';
          ctx.font = '9px Inter, sans-serif';
          ctx.textAlign = 'center';
          ctx.textBaseline = 'middle';
          const label = r.name.length > Math.floor(fW / 7) ? r.name.substring(0, Math.floor(fW / 7)) + '..' : r.name;
          ctx.fillText(label, clampX1 + fW / 2, yOff + featureH / 2);
        }
      });
      currentY += TRACK_H;
    });

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
    const padVal = viewRange * 0.02;
    const viewMin = globalMin - padVal + panOffset;
    const viewMax = viewMin + viewRange + padVal * 2;
    const effRange = viewMax - viewMin;
    return { viewMin, effRange, w: rect.width, h: rect.height };
  };

  const handleWheel = (e: React.WheelEvent) => {
    e.preventDefault();
    const factor = e.deltaY > 0 ? 0.8 : 1.25;
    setZoom(prev => Math.max(0.5, Math.min(10000000, prev * factor)));
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
    const scrollTop = containerRef.current ? containerRef.current.scrollTop : 0;
    const x = e.clientX - rect.left;
    const y = e.clientY - rect.top + scrollTop;

    const bp = p.viewMin + ((x - LABEL_W) / (p.w - LABEL_W)) * p.effRange;
    setCursorBp(bp);

    if (isDragging.current) {
      const dx = e.clientX - dragStartX.current;
      const bpPerPx = p.effRange / (p.w - LABEL_W);
      setPanOffset(panStartVal.current - dx * bpPerPx);
      return;
    }

    const gcOffsetH = gcProfileNames.length * 60;
    const trackIdx = Math.floor((y - RULER_H - gcOffsetH) / TRACK_H);
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
    setZoom(Math.min(10000000, Math.max(0.5, newZoom)));
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
        <button onClick={() => setZoom(prev => Math.min(10000000, prev * 2))}>+ Zoom In</button>
        <button onClick={() => setZoom(prev => Math.max(0.5, prev / 2))}>- Zoom Out</button>
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
              <span className="tt-label">Strand:</span> {formatStrand(hoveredRegion.strand)}<br />
              <span className="tt-label">Length:</span> {(hoveredRegion.end - hoveredRegion.start).toLocaleString()} bp
              {hoveredRegion.motifEvidence && (
                <>
                  <br /><span className="tt-label">Matrix:</span> {hoveredRegion.motifEvidence.matrixId || hoveredRegion.motifEvidence.matrixAlias}
                  <br /><span className="tt-label">Score:</span> {hoveredRegion.motifEvidence.rawScore.toFixed(3)} ({hoveredRegion.motifEvidence.scorePercent.toFixed(1)}%)
                  {hoveredRegion.motifEvidence.statistics && (
                    <><br /><span className="tt-label">Significance:</span> p={hoveredRegion.motifEvidence.statistics.pValue.toExponential(2)}, q={hoveredRegion.motifEvidence.statistics.qValue.toExponential(2)}</>
                  )}
                  {hoveredRegion.motifEvidence.background && (
                    <><br /><span className="tt-label">Background:</span> {hoveredRegion.motifEvidence.background.mode} ({hoveredRegion.motifEvidence.background.source})</>
                  )}
                  {hoveredRegion.motifEvidence.sourceRegion && (
                    <><br /><span className="tt-label">Source:</span> {hoveredRegion.motifEvidence.sourceRegion.name} +{hoveredRegion.motifEvidence.sourceRegion.relativeStart} bp</>
                  )}
                </>
              )}
              {hoveredRegion.spatialRelation && (
                <>
                  <br /><span className="tt-label">{hoveredRegion.spatialRelation.relation}:</span> {hoveredRegion.spatialRelation.reference.name || hoveredRegion.spatialRelation.reference.type}
                  <br /><span className="tt-label">Distance:</span> {formatBp(hoveredRegion.spatialRelation.distance)} / {formatBp(hoveredRegion.spatialRelation.maximumDistance)} max
                  <br /><span className="tt-label">Overlap:</span> {hoveredRegion.spatialRelation.overlaps ? 'yes' : 'no'}
                </>
              )}
              {hoveredRegion.countEvidence && (
                <>
                  <br /><span className="tt-label">Overlap count:</span> {hoveredRegion.countEvidence.count.toLocaleString()}
                  <br /><span className="tt-label">Counted set:</span> {hoveredRegion.countEvidence.countedSet}
                </>
              )}
              {hoveredRegion.consensusEvidence && (
                <>
                  <br /><span className="tt-label">Consensus:</span> {hoveredRegion.consensusEvidence.observedSupport}/{hoveredRegion.consensusEvidence.inputSets.length} sets (minimum {hoveredRegion.consensusEvidence.minimumSupport})
                  {hoveredRegion.consensusEvidence.minimumReciprocalOverlapPercent !== undefined && <><br /><span className="tt-label">Reciprocal overlap:</span> ≥ {hoveredRegion.consensusEvidence.minimumReciprocalOverlapPercent}%</>}
                  <br /><span className="tt-label">Anchor:</span> {hoveredRegion.consensusEvidence.anchorSet}
                </>
              )}
              {hoveredRegion.trackEvidence && (
                <>
                  <br /><span className="tt-label">Track:</span> {hoveredRegion.trackEvidence.trackAlias} ({hoveredRegion.trackEvidence.format})
                  <br /><span className="tt-label">Evidence:</span> {hoveredRegion.trackEvidence.evidenceClass}
                  {hoveredRegion.trackEvidence.assay && <><br /><span className="tt-label">Assay:</span> {hoveredRegion.trackEvidence.assay}</>}
                  {hoveredRegion.trackEvidence.sample && <><br /><span className="tt-label">Sample:</span> {hoveredRegion.trackEvidence.sample}</>}
                  {hoveredRegion.trackEvidence.condition && <><br /><span className="tt-label">Condition:</span> {hoveredRegion.trackEvidence.condition}</>}
                  {hoveredRegion.trackEvidence.replicate && <><br /><span className="tt-label">Replicate:</span> {hoveredRegion.trackEvidence.replicate}</>}
                  {hoveredRegion.trackEvidence.control && <><br /><span className="tt-label">Control:</span> {hoveredRegion.trackEvidence.control}</>}
                  {hoveredRegion.trackEvidence.score !== undefined && <><br /><span className="tt-label">Track score:</span> {hoveredRegion.trackEvidence.score}</>}
                  {hoveredRegion.trackEvidence.signalValue !== undefined && <><br /><span className="tt-label">Signal:</span> {hoveredRegion.trackEvidence.signalValue}</>}
                  {hoveredRegion.trackEvidence.peakPosition !== undefined && <><br /><span className="tt-label">Summit:</span> {hoveredRegion.trackEvidence.peakPosition.toLocaleString()}</>}
                </>
              )}
              {hoveredRegion.overlapEvidence && hoveredRegion.overlapEvidence.length > 0 && (
                <>
                  <br /><span className="tt-label">Overlap support:</span> {hoveredRegion.overlapEvidence.length} reference{hoveredRegion.overlapEvidence.length === 1 ? '' : 's'}
                  <br /><span className="tt-label">Matched:</span> {hoveredRegion.overlapEvidence.map(item => `${item.reference.name || item.reference.type}${item.trackEvidence ? ` [${item.trackEvidence.evidenceClass}${item.trackEvidence.replicate ? `/${item.trackEvidence.replicate}` : ''}]` : ''}${item.supportingEvidence?.length ? ` + ${item.supportingEvidence.length} nested support` : ''}`).join(', ')}
                </>
              )}
              {hoveredRegion.moduleEvidence && (
                <>
                  <br /><span className="tt-label">Module spacing:</span> {formatBp(hoveredRegion.moduleEvidence.spacing.observed)} ({formatBp(hoveredRegion.moduleEvidence.spacing.minimum)}-{formatBp(hoveredRegion.moduleEvidence.spacing.maximum)})
                  <br /><span className="tt-label">Order / orientation:</span> {hoveredRegion.moduleEvidence.order.observed} / {hoveredRegion.moduleEvidence.orientation.observed}
                  <br /><span className="tt-label">Members:</span> {hoveredRegion.moduleEvidence.members.map(member => member.name).join(' + ')}
                </>
              )}
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
                <span className="detail-value">{formatStrand(selectedRegion.strand)}</span>
              </div>
              {selectedRegion.motifEvidence && (
                <>
                  <div className="detail-field">
                    <span className="detail-label">Matrix</span>
                    <span className="detail-value">{selectedRegion.motifEvidence.matrixId || selectedRegion.motifEvidence.matrixAlias} ({selectedRegion.motifEvidence.matrixName})</span>
                  </div>
                  <div className="detail-field">
                    <span className="detail-label">PWM score</span>
                    <span className="detail-value">{selectedRegion.motifEvidence.rawScore.toFixed(3)} / {selectedRegion.motifEvidence.scorePercent.toFixed(1)}%</span>
                  </div>
                  {selectedRegion.motifEvidence.statistics && (
                    <div className="detail-field">
                      <span className="detail-label">Statistical evidence</span>
                      <span className="detail-value">
                        p {selectedRegion.motifEvidence.statistics.pValue.toExponential(4)} · q {selectedRegion.motifEvidence.statistics.qValue.toExponential(4)}<br />
                        {selectedRegion.motifEvidence.statistics.testedPositions.toLocaleString()} tested position-strands · {selectedRegion.motifEvidence.statistics.multipleTestingMethod}
                      </span>
                    </div>
                  )}
                  {selectedRegion.motifEvidence.background && (
                    <div className="detail-field">
                      <span className="detail-label">Background</span>
                      <span className="detail-value">
                        {selectedRegion.motifEvidence.background.mode} from {selectedRegion.motifEvidence.background.source} ({selectedRegion.motifEvidence.background.strandPolicy})<br />
                        A {selectedRegion.motifEvidence.background.A.toFixed(3)} · C {selectedRegion.motifEvidence.background.C.toFixed(3)} · G {selectedRegion.motifEvidence.background.G.toFixed(3)} · T {selectedRegion.motifEvidence.background.T.toFixed(3)}
                      </span>
                    </div>
                  )}
                  {selectedRegion.motifEvidence.sourceRegion && (
                    <div className="detail-field">
                      <span className="detail-label">Source region</span>
                      <span className="detail-value">{selectedRegion.motifEvidence.sourceRegion.name}, +{selectedRegion.motifEvidence.sourceRegion.relativeStart} bp</span>
                    </div>
                  )}
                </>
              )}
              {selectedRegion.spatialRelation && (
                <>
                  <div className="detail-field">
                    <span className="detail-label">Spatial relation</span>
                    <span className="detail-value">
                      {selectedRegion.spatialRelation.relation} {selectedRegion.spatialRelation.referenceSet}<br />
                      {selectedRegion.spatialRelation.reference.name || selectedRegion.spatialRelation.reference.type} · {selectedRegion.spatialRelation.reference.chr}:{selectedRegion.spatialRelation.reference.start.toLocaleString()}-{selectedRegion.spatialRelation.reference.end.toLocaleString()}
                    </span>
                  </div>
                  <div className="detail-field">
                    <span className="detail-label">Interval distance</span>
                    <span className="detail-value">
                      {formatBp(selectedRegion.spatialRelation.distance)} (limit {formatBp(selectedRegion.spatialRelation.maximumDistance)}) · {selectedRegion.spatialRelation.overlaps ? 'overlapping' : 'not overlapping'}
                    </span>
                  </div>
                </>
              )}
              {selectedRegion.countEvidence && (
                <div className="detail-field">
                  <span className="detail-label">Overlap count</span>
                  <span className="detail-value">
                    {selectedRegion.countEvidence.count.toLocaleString()} from {selectedRegion.countEvidence.countedSet}<br />
                    {selectedRegion.countEvidence.relation} in {selectedRegion.countEvidence.containerSet}
                  </span>
                </div>
              )}
              {selectedRegion.trackEvidence && (
                <>
                  <div className="detail-field">
                    <span className="detail-label">Imported track</span>
                    <span className="detail-value">
                      {selectedRegion.trackEvidence.trackAlias} · {selectedRegion.trackEvidence.format}<br />
                      {selectedRegion.trackEvidence.source}<br />
                      {selectedRegion.trackEvidence.evidenceClass}
                      {selectedRegion.trackEvidence.assay ? ` · ${selectedRegion.trackEvidence.assay}` : ''}
                      {selectedRegion.trackEvidence.sample ? ` · ${selectedRegion.trackEvidence.sample}` : ''}
                      {selectedRegion.trackEvidence.condition ? ` · condition ${selectedRegion.trackEvidence.condition}` : ''}
                      {selectedRegion.trackEvidence.replicate ? ` · replicate ${selectedRegion.trackEvidence.replicate}` : ''}
                      {selectedRegion.trackEvidence.control ? ` · control ${selectedRegion.trackEvidence.control}` : ''}
                    </span>
                  </div>
                  <div className="detail-field">
                    <span className="detail-label">Track evidence</span>
                    <span className="detail-value">
                      score {selectedRegion.trackEvidence.score ?? 'N/A'}
                      {selectedRegion.trackEvidence.signalValue !== undefined ? ` · signal ${selectedRegion.trackEvidence.signalValue}` : ''}
                      {selectedRegion.trackEvidence.minusLog10PValue !== undefined ? ` · -log10(p) ${selectedRegion.trackEvidence.minusLog10PValue}` : ''}
                      {selectedRegion.trackEvidence.minusLog10QValue !== undefined ? ` · -log10(q) ${selectedRegion.trackEvidence.minusLog10QValue}` : ''}
                      {selectedRegion.trackEvidence.peakPosition !== undefined ? ` · summit ${selectedRegion.trackEvidence.peakPosition.toLocaleString()}` : ''}
                    </span>
                  </div>
                </>
              )}
              {selectedRegion.consensusEvidence && (
                <div className="detail-field">
                  <span className="detail-label">Coordinate consensus</span>
                  <span className="detail-value">
                    {selectedRegion.consensusEvidence.observedSupport}/{selectedRegion.consensusEvidence.inputSets.length} supporting sets · minimum {selectedRegion.consensusEvidence.minimumSupport}<br />
                    anchor {selectedRegion.consensusEvidence.anchorSet}<br />
                    {selectedRegion.consensusEvidence.minimumReciprocalOverlapPercent !== undefined && <>reciprocal overlap ≥ {selectedRegion.consensusEvidence.minimumReciprocalOverlapPercent}%<br /></>}
                    {selectedRegion.consensusEvidence.inputSets.join(', ')}
                  </span>
                </div>
              )}
              {selectedRegion.overlapEvidence && selectedRegion.overlapEvidence.length > 0 && (
                <div className="detail-field">
                  <span className="detail-label">Overlap evidence</span>
                  <span className="detail-value">
                    {selectedRegion.overlapEvidence.map((item, index) => (
                      <span key={`${item.referenceSet}-${item.reference.chr}-${item.reference.start}-${index}`} style={{ display: 'block' }}>
                        {item.referenceSet}: {item.reference.name || item.reference.type} · {item.reference.chr}:{item.reference.start.toLocaleString()}-{item.reference.end.toLocaleString()}
                        {item.trackEvidence ? ` · ${item.trackEvidence.evidenceClass}${item.trackEvidence.assay ? ` (${item.trackEvidence.assay})` : ''}` : ''}
                        {item.trackEvidence?.condition ? ` · condition ${item.trackEvidence.condition}` : ''}
                        {item.trackEvidence?.replicate ? ` · replicate ${item.trackEvidence.replicate}` : ''}
                        {item.trackEvidence?.control ? ` · control ${item.trackEvidence.control}` : ''}
                        {item.consensusEvidence ? ` · consensus ${item.consensusEvidence.observedSupport}/${item.consensusEvidence.inputSets.length} (minimum ${item.consensusEvidence.minimumSupport}${item.consensusEvidence.minimumReciprocalOverlapPercent !== undefined ? `; reciprocal overlap ≥ ${item.consensusEvidence.minimumReciprocalOverlapPercent}%` : ''})` : ''}
                        {item.supportingEvidence?.map((support, supportIndex) => (
                          <span key={`${support.referenceSet}-${support.reference.chr}-${support.reference.start}-${supportIndex}`} style={{ display: 'block', paddingLeft: 12 }}>
                            ↳ {support.referenceSet}: {support.reference.name || support.reference.type}
                            {support.trackEvidence ? ` · ${support.trackEvidence.evidenceClass}` : ''}
                            {support.trackEvidence?.condition ? ` · condition ${support.trackEvidence.condition}` : ''}
                            {support.trackEvidence?.replicate ? ` · replicate ${support.trackEvidence.replicate}` : ''}
                            {support.trackEvidence?.control ? ` · control ${support.trackEvidence.control}` : ''}
                          </span>
                        ))}
                      </span>
                    ))}
                  </span>
                </div>
              )}
              {selectedRegion.moduleEvidence && (
                <>
                  <div className="detail-field">
                    <span className="detail-label">Module constraints</span>
                    <span className="detail-value">
                      {formatBp(selectedRegion.moduleEvidence.spacing.observed)} spacing (allowed {formatBp(selectedRegion.moduleEvidence.spacing.minimum)}-{formatBp(selectedRegion.moduleEvidence.spacing.maximum)})<br />
                      order {selectedRegion.moduleEvidence.order.observed} ({selectedRegion.moduleEvidence.order.policy}) · orientation {selectedRegion.moduleEvidence.orientation.observed} ({selectedRegion.moduleEvidence.orientation.policy})
                    </span>
                  </div>
                  <div className="detail-field">
                    <span className="detail-label">Module members</span>
                    <span className="detail-value">
                      {selectedRegion.moduleEvidence.members.map((member, index) => (
                        <React.Fragment key={`${member.sourceSet}-${member.chr}-${member.start}-${member.end}-${index}`}>
                          {index > 0 && <br />}
                          {member.sourceSet}: {member.name} · {member.chr}:{member.start.toLocaleString()}-{member.end.toLocaleString()} · {member.strand}{member.motifEvidence ? ` · ${member.motifEvidence.matrixId || member.motifEvidence.matrixAlias}` : ''}
                        </React.Fragment>
                      ))}
                    </span>
                  </div>
                </>
              )}
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
