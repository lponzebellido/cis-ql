import React, { useRef, useEffect, useState } from 'react';

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
}

const COLORS = [
  '#4ade80', 
  '#60a5fa', 
  '#f472b6', 
  '#fbbf24', 
  '#c084fc', 
  '#2dd4bf', 
];

export const TrackViewer: React.FC<TrackViewerProps> = ({ results }) => {
  const canvasRef = useRef<HTMLCanvasElement>(null);
  const containerRef = useRef<HTMLDivElement>(null);
  const [hoveredRegion, setHoveredRegion] = useState<GenomicRegion | null>(null);
  const [mousePos, setMousePos] = useState({ x: 0, y: 0 });

  
  let minPos = Infinity;
  let maxPos = 0;

  const trackNames = Object.keys(results);
  const allRegions: (GenomicRegion & { trackIdx: number })[] = [];

  trackNames.forEach((trackName, idx) => {
    results[trackName].forEach(r => {
      if (r.start < minPos) minPos = r.start;
      if (r.end > maxPos) maxPos = r.end;
      allRegions.push({ ...r, trackIdx: idx });
    });
  });

  if (minPos === Infinity) {
    minPos = 0;
    maxPos = 1000; 
  }

  
  const range = maxPos - minPos || 1000;
  const padding = range * 0.05;
  const drawMin = Math.max(0, minPos - padding);
  const drawMax = maxPos + padding;
  const drawRange = drawMax - drawMin;

  useEffect(() => {
    const canvas = canvasRef.current;
    const container = containerRef.current;
    if (!canvas || !container) return;

    const ctx = canvas.getContext('2d');
    if (!ctx) return;

    
    const rect = container.getBoundingClientRect();
    canvas.width = rect.width * window.devicePixelRatio;
    canvas.height = rect.height * window.devicePixelRatio;
    ctx.scale(window.devicePixelRatio, window.devicePixelRatio);
    canvas.style.width = `${rect.width}px`;
    canvas.style.height = `${rect.height}px`;

    const width = rect.width;
    const height = rect.height;

    
    ctx.fillStyle = '#1e1e2e'; 
    ctx.fillRect(0, 0, width, height);

    if (trackNames.length === 0) {
      ctx.fillStyle = '#a6adc8';
      ctx.font = '14px Inter, sans-serif';
      ctx.textAlign = 'center';
      ctx.fillText('No regions to display', width / 2, height / 2);
      return;
    }

    const trackHeight = height / (trackNames.length + 1); 

    
    ctx.strokeStyle = '#45475a';
    ctx.lineWidth = 1;
    ctx.beginPath();
    const axisY = height - 30;
    ctx.moveTo(0, axisY);
    ctx.lineTo(width, axisY);
    ctx.stroke();

    ctx.fillStyle = '#a6adc8';
    ctx.font = '12px Inter, sans-serif';
    ctx.textAlign = 'center';
    
    
    for (let i = 0; i <= 5; i++) {
      const x = (i / 5) * width;
      const bp = Math.floor(drawMin + (i / 5) * drawRange);
      ctx.beginPath();
      ctx.moveTo(x, axisY - 5);
      ctx.lineTo(x, axisY + 5);
      ctx.stroke();
      ctx.fillText(`${bp.toLocaleString()} bp`, x, axisY + 20);
    }

    
    trackNames.forEach((trackName, idx) => {
      const yCenter = (idx + 0.5) * trackHeight;
      
      
      ctx.fillStyle = '#cdd6f4';
      ctx.textAlign = 'left';
      ctx.fillText(trackName, 10, yCenter - 15);
      
      
      ctx.strokeStyle = '#313244';
      ctx.beginPath();
      ctx.moveTo(0, yCenter);
      ctx.lineTo(width, yCenter);
      ctx.stroke();

      const color = COLORS[idx % COLORS.length];
      
      
      results[trackName].forEach(r => {
        const xStart = ((r.start - drawMin) / drawRange) * width;
        const xEnd = ((r.end - drawMin) / drawRange) * width;
        const featWidth = Math.max(1, xEnd - xStart); 
        
        ctx.fillStyle = color;
        
        const boxY = r.strand === '+' ? yCenter - 12 : yCenter + 2;
        ctx.fillRect(xStart, boxY, featWidth, 10);
      });
    });

  }, [results, drawMin, drawMax, drawRange, trackNames]);

  const handleMouseMove = (e: React.MouseEvent<HTMLCanvasElement>) => {
    const rect = canvasRef.current?.getBoundingClientRect();
    if (!rect || trackNames.length === 0) return;

    const x = e.clientX - rect.left;
    const y = e.clientY - rect.top;
    
    
    const bp = drawMin + (x / rect.width) * drawRange;
    
    
    const trackHeight = rect.height / (trackNames.length + 1);
    const trackIdx = Math.floor(y / trackHeight);

    let found: GenomicRegion | null = null;
    
    if (trackIdx >= 0 && trackIdx < trackNames.length) {
      const trackName = trackNames[trackIdx];
      
      
      const paddingBp = drawRange * 0.005; 
      for (const r of results[trackName]) {
        if (bp >= r.start - paddingBp && bp <= r.end + paddingBp) {
          found = r;
          break;
        }
      }
    }

    setHoveredRegion(found);
    setMousePos({ x: e.clientX, y: e.clientY });
  };

  const handleMouseLeave = () => {
    setHoveredRegion(null);
  };

  return (
    <div style={{ position: 'relative', width: '100%', height: '100%' }} ref={containerRef}>
      <canvas 
        ref={canvasRef} 
        onMouseMove={handleMouseMove}
        onMouseLeave={handleMouseLeave}
        style={{ cursor: hoveredRegion ? 'pointer' : 'default' }}
      />
      
      {hoveredRegion && (
        <div style={{
          position: 'fixed',
          left: mousePos.x + 15,
          top: mousePos.y + 15,
          backgroundColor: 'rgba(30, 30, 46, 0.95)',
          border: '1px solid #45475a',
          borderRadius: '6px',
          padding: '12px',
          color: '#cdd6f4',
          fontFamily: 'Inter, sans-serif',
          fontSize: '12px',
          pointerEvents: 'none',
          boxShadow: '0 4px 6px rgba(0,0,0,0.3)',
          zIndex: 1000
        }}>
          <strong style={{ color: '#89b4fa', fontSize: '14px' }}>{hoveredRegion.name}</strong><br/>
          <div style={{ marginTop: '4px' }}>
            <span style={{ color: '#a6adc8' }}>Type:</span> {hoveredRegion.type}<br/>
            <span style={{ color: '#a6adc8' }}>Location:</span> {hoveredRegion.chr}:{hoveredRegion.start}-{hoveredRegion.end} ({hoveredRegion.strand})<br/>
            {hoveredRegion.sequence && (
              <div style={{ marginTop: '6px', fontFamily: 'monospace', color: '#a6e3a1' }}>
                {hoveredRegion.sequence}
              </div>
            )}
          </div>
        </div>
      )}
    </div>
  );
};
