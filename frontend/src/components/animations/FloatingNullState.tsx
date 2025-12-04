import React, { useRef, useState, useEffect, useMemo } from 'react';
import { motion, AnimatePresence } from 'framer-motion';

// Simplified: static tick ring with spiky waves at random angles.
// Every ~3–4s, spawn four provider items at random angles (>=60° apart)
// and sizes; the angles drive a smooth wave on ticks.
type FloatingNullStateTiming = {
  pushDurationMs?: number;
  relaxDurationMs?: number;
  itemHoldMs?: number;
  idleMinMs?: number;
  idleJitterMs?: number;
  staggerStepMs?: number;
  staggerJitterMs?: number;
};

const FloatingNullStateComponent: React.FC<{ hideIcons?: boolean; timing?: FloatingNullStateTiming; compact?: boolean }> = ({ hideIcons = false, timing, compact = false }) => {
  const tickCount = 96;
  const innerRadius = 36; // ring inner radius (in SVG viewBox units)
  const outerRadius = 44; // ring outer radius

  // Smooth wave parameters (deg-based)
  const WAVE_WIDTH_DEG = 15;         // ~8 ticks wide at 96 ticks for smooth gradient
  const WAVE_EXP = 2.5;              // gentler exponent => smoother bell curve
  const LENGTH_BOOST = 7;            // extra outer radius when wave is strongest
  const OPACITY_BASE = 0.20;         // baseline stroke opacity
  const OPACITY_PEAK = 1.0;          // peak opacity at wave center
  const TICK_TRANSLATE_UNITS = 1.0;  // outward translation in viewBox units at peak
  const TICK_ANGLE_DEG = 360 / 96;   // 3.75° per tick
  const ITEM_PUSH_MIN_PX = 12;        // broaden variation slightly
  const ITEM_PUSH_MAX_PX = 18;
  const DEBUG = false;               // toggle to show a red dot at exact polar point
  const FADE_OUT_START_AMP = 0.45;  // begin fading before items reach base radius
  const BASE_CENTER_OFFSET_PX = 6;  // constant outward offset so items never drift onto ticks
  const EXCLUDED_ARCS: Array<{ center: number; halfWidth: number }> = [
    { center: 0, halfWidth: 30 },
    { center: 180, halfWidth: 30 },
  ];

  // Timing (overridable via props.timing)
  const {
    pushDurationMs,
    relaxDurationMs,
    itemHoldMs,
    idleMinMs,
    idleJitterMs,
    staggerStepMs,
    staggerJitterMs,
  } = {
    pushDurationMs: 1200,
    relaxDurationMs: 700,
    itemHoldMs: 300,
    idleMinMs: 1000,
    idleJitterMs: 500,
    staggerStepMs: 150,
    staggerJitterMs: 30,
    ...(timing || {}),
  };

  const PUSH_DURATION_MS = pushDurationMs;
  const RELAX_DURATION_MS = relaxDurationMs;
  const ITEM_HOLD_MS = itemHoldMs;          // hold at peak before fading items out
  const IDLE_MIN_MS = idleMinMs;            // rest between spikes
  const IDLE_JITTER_MS = idleJitterMs;      // + up to 1s
  const STAGGER_STEP_MS = staggerStepMs;    // base delta-T between items in same batch
  const STAGGER_JITTER_MS = staggerJitterMs; // +/- jitter per item

  // Helper: shortest angular distance accounting for 0/360 wrap
  function angleDiff(a: number, b: number): number {
    let d = a - b;
    while (d > 180) d -= 360;
    while (d < -180) d += 360;
    return d;
  }

  // Provider icon catalog
  const providers: Array<{ src: string; alt: string; offsetPx?: { dx: number; dy: number }; radialOffsetPx?: number }> = [
    { src: '/images/ChemOne.png', alt: 'Chemistry 1' },
    { src: '/images/ChemTwo.png', alt: 'Chemistry 2' },
    { src: '/images/ChemThree.png', alt: 'Chemistry 3' },
    { src: '/images/ChemFour.png', alt: 'Chemistry 4' },
    { src: '/images/ChemFive.png', alt: 'Chemistry 5' },
    { src: '/images/ChemSix.png', alt: 'Chemistry 6' },
    { src: '/images/ChemSeven.png', alt: 'Chemistry 7' },
    { src: '/images/ChemEight.png', alt: 'Chemistry 8' },
    { src: '/images/ChemNine.png', alt: 'Chemistry 9' },
  ];

  // Convert angle+radius to CSS percentage positions relative to 100x100 viewbox
  // angleDeg: 0 = top, clockwise (same as tick/wave coordinate system)
  function posFromAngle(angleDeg: number, radius: number) {
    const angle = (angleDeg * Math.PI) / 180;
    const cx = 50; const cy = 50;
    const angleOffset = angle - Math.PI / 2; // offset so 0° = top
    const x = cx + radius * Math.cos(angleOffset);
    const y = cy + radius * Math.sin(angleOffset);
    return { left: `${x}%`, top: `${y}%` };
  }

  // Active wave peaks (with per-wave delays) and active items
  const [activeWaves, setActiveWaves] = useState<Array<{ angleDeg: number; delayMs: number }>>([]);
  const [activeItems, setActiveItems] = useState<Array<{ id: number; src: string; alt: string; size: number; angleDeg: number; pushPx: number; bufferPx: number; radiusVB: number; delayMs: number; offsetPx?: { dx: number; dy: number }; radialOffsetPx?: number }>>([]);
  const [nowMs, setNowMs] = useState(0);
  const nextIdRef = useRef(1);
  const animRafRef = useRef<number | null>(null);
  const timeoutsRef = useRef<number[]>([]);
  const stageRef = useRef<HTMLDivElement | null>(null);
  const [pxToVb, setPxToVb] = useState<number>(100 / 300); // default assuming 300px stage
  const batchStartRef = useRef<number>(0);
  // Viewport-height driven scaling and spacing
  const [viewportH, setViewportH] = useState<number>(typeof window !== 'undefined' ? window.innerHeight : 800);
  useEffect(() => {
    const onResize = () => { try { setViewportH(window.innerHeight); } catch (_) {} };
    window.addEventListener('resize', onResize);
    return () => window.removeEventListener('resize', onResize);
  }, []);
  const SCALE_MIN_VH = 600;
  const SCALE_MAX_VH = 1100;
  const vhT = Math.max(0, Math.min(1, (viewportH - SCALE_MIN_VH) / (SCALE_MAX_VH - SCALE_MIN_VH)));
  const scaleFactor = 1 + 0.5 * vhT; // 1.0 .. 1.5
  const marginTopBase = Math.round(8 + 12 * vhT);  // 8 .. 20
  const marginBottomBase = Math.round(16 + 24 * vhT); // 16 .. 40
  const marginTopPx = compact ? Math.max(4, Math.round(marginTopBase * 0.6)) : marginTopBase;
  const marginBottomPx = compact ? Math.max(8, Math.round(marginBottomBase * 0.6)) : marginBottomBase;

  // keep px->viewBox conversion accurate to stage size and applied scale
  useEffect(() => {
    const update = () => {
      if (!stageRef.current) return;
      const rect = stageRef.current.getBoundingClientRect();
      if (rect.width > 0) setPxToVb(100 / rect.width);
    };
    update();
    let ro: ResizeObserver | null = null;
    if (typeof ResizeObserver !== 'undefined') {
      ro = new ResizeObserver(update);
      if (stageRef.current) ro.observe(stageRef.current);
    } else {
      window.addEventListener('resize', update);
    }
    return () => {
      if (ro) ro.disconnect();
      else window.removeEventListener('resize', update);
    };
  }, [scaleFactor]);

  function randomAngle(): number {
    return Math.random() * 360;
  }

  function snapToNearestTick(angleDeg: number): number {
    // Snap to nearest tick angle for perfect alignment
    const tickIndex = Math.round(angleDeg / TICK_ANGLE_DEG);
    return (tickIndex % tickCount) * TICK_ANGLE_DEG;
  }

  function angleDistanceDeg(a: number, b: number): number {
    return Math.abs(angleDiff(a, b));
  }

  function isAngleExcluded(angleDeg: number): boolean {
    for (const arc of EXCLUDED_ARCS) {
      if (angleDistanceDeg(angleDeg, arc.center) < arc.halfWidth) return true;
    }
    return false;
  }

  function sampleAngles(count: number, minSeparationDeg = 60): number[] {
    // Sample N random angles then snap to tick positions with minimum separation
    const angles: number[] = [];
    const maxAttempts = 100;
    
    for (let n = 0; n < count; n++) {
      let found = false;
      for (let attempt = 0; attempt < maxAttempts; attempt++) {
        const candidate = snapToNearestTick(randomAngle());
        // Skip excluded arcs, then check separation from existing angles
        const valid = !isAngleExcluded(candidate) &&
          angles.every(existing => angleDistanceDeg(candidate, existing) >= minSeparationDeg);
        if (valid) {
          angles.push(candidate);
          found = true;
          break;
        }
      }
      // Fallback: distribute evenly around the circle with small jitter
      if (!found) {
        const baseAngle = (n * 360) / count;
        const jitter = (Math.random() * 20) - 10; // [-10, +10]
        let fb = snapToNearestTick((baseAngle + jitter + 360) % 360);
        // bump by one tick step until we leave excluded arcs, with safety bound
        let safety = 0;
        while (isAngleExcluded(fb) && safety < tickCount) {
          fb = snapToNearestTick((fb + TICK_ANGLE_DEG + 360) % 360);
          safety++;
        }
        angles.push(fb);
      }
    }
    return angles;
  }

  function randomSize(min = 16, max = 42): number {
    return Math.round(min + Math.random() * (max - min));
  }

  function sampleProviders(count: number): typeof providers[number][] {
    if (providers.length === 0) return [];
    // Shuffle and pick first N (with replacement if needed)
    const shuffled = [...providers].sort(() => Math.random() - 0.5);
    const result: typeof providers[number][] = [];
    for (let i = 0; i < count; i++) {
      result.push(shuffled[i % shuffled.length]);
    }
    return result;
  }

  function spawnCycle() {
    const itemCount = 4;
    const angles = sampleAngles(itemCount, 60);
    const providerList = sampleProviders(itemCount);
    
    // Per-item/wave delays to create organic staggering (controllable)
    const delays = angles.map((_, i) => {
      const jitter = Math.round((Math.random() - 0.5) * 2 * STAGGER_JITTER_MS);
      return Math.max(0, i * STAGGER_STEP_MS + jitter);
    });

    const items = angles.map((angleDeg, idx) => {
      const size = randomSize(12, 28);
      const bufferPx = 1 + Math.floor(Math.random() * 3); // 1..3 px
      const pushPx = ITEM_PUSH_MIN_PX + Math.random() * (ITEM_PUSH_MAX_PX - ITEM_PUSH_MIN_PX);
      
      // Polar coordinate positioning:
      // Start position: outerRadius (tick outer edge) + bufferPx + itemRadius (to align item inner edge with tick outer)
      // itemRadius in viewBox units = (size + 15) / 2 converted to VB
      const itemRadiusPx = (size + 15) / 2;
      const itemRadiusVB = itemRadiusPx * pxToVb;
      const radiusVB = outerRadius + bufferPx * pxToVb + itemRadiusVB;
      
      return {
        id: nextIdRef.current++,
        src: providerList[idx].src,
        alt: providerList[idx].alt,
        size,
        angleDeg,
        pushPx,
        bufferPx,
        radiusVB,
        delayMs: delays[idx],
        offsetPx: providerList[idx].offsetPx,
        radialOffsetPx: providerList[idx].radialOffsetPx,
      };
    });
    setActiveWaves(angles.map((angleDeg, i) => ({ angleDeg, delayMs: delays[i] })));
    setActiveItems(items);

    // Start batch clock
    batchStartRef.current = performance.now();
    setNowMs(batchStartRef.current);
    const raf = () => {
      setNowMs(performance.now());
      animRafRef.current = requestAnimationFrame(raf);
    };
    if (animRafRef.current) cancelAnimationFrame(animRafRef.current);
    animRafRef.current = requestAnimationFrame(raf);

    // Schedule end-of-cycle and next spawn after the last delayed wave completes
    const maxDelay = Math.max(...delays, 0);
    const cycleMs = PUSH_DURATION_MS + ITEM_HOLD_MS + RELAX_DURATION_MS;
    const totalMs = maxDelay + cycleMs;
    const endTimer = window.setTimeout(() => {
      // stop raf, clear items/waves, then idle and respawn
      if (animRafRef.current) cancelAnimationFrame(animRafRef.current);
      animRafRef.current = null;
      setActiveItems([]);
      setActiveWaves([]);
      const idleDelay = IDLE_MIN_MS + Math.random() * IDLE_JITTER_MS;
      const nextTimer = window.setTimeout(spawnCycle, idleDelay);
      timeoutsRef.current.push(nextTimer);
    }, totalMs);
    timeoutsRef.current.push(endTimer);
  }

  useEffect(() => {
    spawnCycle();
    return () => {
      if (animRafRef.current) cancelAnimationFrame(animRafRef.current);
      for (const t of timeoutsRef.current) window.clearTimeout(t);
      timeoutsRef.current = [];
    };
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, []);

  function easeInOutCubic(t: number) {
    return t < 0.5 ? 4 * t * t * t : 1 - Math.pow(-2 * t + 2, 3) / 2;
  }

  // Interpolate between two colors (hex format) based on a value t (0-1)
  function interpolateColor(color1: string, color2: string, t: number): string {
    // Parse hex colors
    const c1 = {
      r: parseInt(color1.slice(1, 3), 16),
      g: parseInt(color1.slice(3, 5), 16),
      b: parseInt(color1.slice(5, 7), 16)
    };
    const c2 = {
      r: parseInt(color2.slice(1, 3), 16),
      g: parseInt(color2.slice(3, 5), 16),
      b: parseInt(color2.slice(5, 7), 16)
    };
    // Interpolate
    const r = Math.round(c1.r + (c2.r - c1.r) * t);
    const g = Math.round(c1.g + (c2.g - c1.g) * t);
    const b = Math.round(c1.b + (c2.b - c1.b) * t);
    return `#${r.toString(16).padStart(2, '0')}${g.toString(16).padStart(2, '0')}${b.toString(16).padStart(2, '0')}`;
  }

  // Local per-wave phase (0..1 during push/hold/relax, else 0)
  function localPhaseFor(nowAbs: number, delayMs: number): number {
    const t = nowAbs - (batchStartRef.current + delayMs);
    if (t <= 0) return 0;
    if (t <= PUSH_DURATION_MS) return easeInOutCubic(t / PUSH_DURATION_MS);
    if (t <= PUSH_DURATION_MS + ITEM_HOLD_MS) return 1;
    if (t <= PUSH_DURATION_MS + ITEM_HOLD_MS + RELAX_DURATION_MS) {
      const rt = (t - PUSH_DURATION_MS - ITEM_HOLD_MS) / RELAX_DURATION_MS;
      return easeInOutCubic(1 - rt);
    }
    return 0;
  }

  // Outward + inward phase for orbit items:
  // - Push: 0 -> 1
  // - Hold: 1
  // - Relax: 1 -> 0.5 (collapse halfway back inwards)
  function outwardInwardPhaseFor(nowAbs: number, delayMs: number): number {
    const t = nowAbs - (batchStartRef.current + delayMs);
    if (t <= 0) return 0;
    if (t <= PUSH_DURATION_MS) return easeInOutCubic(t / PUSH_DURATION_MS);
    if (t <= PUSH_DURATION_MS + ITEM_HOLD_MS) return 1;
    if (t <= PUSH_DURATION_MS + ITEM_HOLD_MS + RELAX_DURATION_MS) {
      const rt = (t - PUSH_DURATION_MS - ITEM_HOLD_MS) / RELAX_DURATION_MS; // 0..1
      const eased = easeInOutCubic(rt);
      return 1 - 0.5 * eased; // collapse to 0.5 at end of relax
    }
    return 0.5;
  }

  // 0..1 progress strictly across relax window; 0 outside relax
  function relaxProgressFor(nowAbs: number, delayMs: number): number {
    const t = nowAbs - (batchStartRef.current + delayMs);
    const start = PUSH_DURATION_MS + ITEM_HOLD_MS;
    if (t <= start) return 0;
    if (t >= start + RELAX_DURATION_MS) return 1;
    return (t - start) / RELAX_DURATION_MS;
  }

  // Returns 0..1 wave strength at a given angle (deg) and absolute time,
  // matching the exact logic used by the SVG tick rendering.
  function waveStrengthAtAngleDeg(angleDeg: number, nowAbs: number): number {
    let strength = 0;
    for (const wave of activeWaves) {
      const d = Math.abs(angleDiff(angleDeg, wave.angleDeg));
      const normalized = d / WAVE_WIDTH_DEG; // 0 at center, 1 at edge
      const falloff = Math.max(0, 1 - normalized);
      const amp = localPhaseFor(nowAbs, wave.delayMs); // 0..1 easing
      const contrib = Math.pow(falloff, WAVE_EXP) * amp;
      if (contrib > strength) strength = contrib;
    }
    return strength;
  }

  // Memoize tick calculations to avoid recalculating on every render
  const tickElements = useMemo(() => {
    const BASE_COLOR = '#6B7C78';   // Darker gray-blue for surrounding ticks
    // Biomarker purple on white brand color for spiking ticks
    const PEAK_COLOR = '#8C6FC1';   // matches bg-biomarker-purple-on-white
    // Global glow state: whenever any wave is active (push/hold/relax), glow the whole ring
    const glowOn = activeWaves.some(wave => localPhaseFor(nowMs, wave.delayMs) > 0);
    
    return Array.from({ length: tickCount }).map((_, i) => {
      // Tick angle in degrees (0 = top)
      const tickDeg = (i / tickCount) * 360;
      const angle = (tickDeg * Math.PI) / 180; // convert to radians for rendering
      const cx = 50;
      const cy = 50;
      const angleOffset = angle - Math.PI / 2; // offset so 0° is top
      // Wave influence from current peaks (smooth bell curve)
      const strength = waveStrengthAtAngleDeg(tickDeg, nowMs);
      const extra = LENGTH_BOOST * strength; // length modulation
      const or = outerRadius + extra;
      const opacity = OPACITY_BASE + (OPACITY_PEAK - OPACITY_BASE) * strength;
      // Interpolate color from base to peak based on local wave strength
      const color = interpolateColor(BASE_COLOR, PEAK_COLOR, strength);
      // Gentle outward translation of the whole tick segment
      const translateUnits = TICK_TRANSLATE_UNITS * strength;
      const x1t = cx + (innerRadius + translateUnits) * Math.cos(angleOffset);
      const y1t = cy + (innerRadius + translateUnits) * Math.sin(angleOffset);
      const x2t = cx + (or + translateUnits) * Math.cos(angleOffset);
      const y2t = cy + (or + translateUnits) * Math.sin(angleOffset);
      return (
        <line
          key={i}
          x1={x1t}
          y1={y1t}
          x2={x2t}
          y2={y2t}
          stroke={color}
          strokeWidth="0.35"
          strokeLinecap="butt"
          opacity={opacity}
          // Subtle biomarker-purple glow whenever the ring is in an expansion cycle
          filter={glowOn ? 'url(#tick-glow)' : undefined}
        />
      );
    });
  }, [activeWaves, nowMs, tickCount, innerRadius, outerRadius, WAVE_WIDTH_DEG, WAVE_EXP, LENGTH_BOOST, OPACITY_BASE, OPACITY_PEAK]);

  return (
    <div 
      className="relative w-full flex items-center justify-center" 
      aria-hidden
      style={{ 
        pointerEvents: 'none',
        userSelect: 'none',
        touchAction: 'none',
        marginTop: marginTopPx,
        marginBottom: marginBottomPx
      }}
    >
      {/* Square stage, responsive size */}
      <div
        className="relative"
        style={{ 
          width: 'clamp(120px, 34vw, 240px)', 
          height: 'clamp(120px, 34vw, 240px)', 
          pointerEvents: 'none',
          transform: `scale(${scaleFactor})`,
          transformOrigin: 'center center'
        }}
        ref={stageRef}
      >
        {/* Tick ring */}
        <svg
          className="absolute inset-0 m-auto"
          width="100%"
          height="100%"
          viewBox="0 0 100 100"
          preserveAspectRatio="xMidYMid meet"
          style={{ pointerEvents: 'none' }}
        >
          {/* Brand-aligned glow used for spiking ticks */}
          <defs>
            <filter id="tick-glow" x="-50%" y="-50%" width="200%" height="200%">
              {/* Blur the tick strokes to get a soft mask */}
              <feGaussianBlur in="SourceGraphic" stdDeviation="0.8" result="blur" />
              {/* Color the blur with biomarker purple-on-white */}
              <feFlood floodColor="#8C6FC1" floodOpacity="0.45" result="purple" />
              <feComposite in="purple" in2="blur" operator="in" result="purpleGlow" />
              <feMerge>
                {/* Glow underneath */}
                <feMergeNode in="purpleGlow" />
                {/* Original ticks on top */}
                <feMergeNode in="SourceGraphic" />
              </feMerge>
            </filter>
          </defs>
          {tickElements}
          {DEBUG && (
            <g>
              {EXCLUDED_ARCS.map((arc, idx) => {
                const cx = 50; const cy = 50;
                const r = (innerRadius + outerRadius) / 2;
                const startDeg = (arc.center - arc.halfWidth + 360) % 360;
                const endDeg = (arc.center + arc.halfWidth + 360) % 360;
                const startRad = (startDeg * Math.PI) / 180;
                const endRad = (endDeg * Math.PI) / 180;
                const startAngle = startRad - Math.PI / 2;
                const endAngle = endRad - Math.PI / 2;
                const x1 = cx + r * Math.cos(startAngle);
                const y1 = cy + r * Math.sin(startAngle);
                const x2 = cx + r * Math.cos(endAngle);
                const y2 = cy + r * Math.sin(endAngle);
                const d = `M ${x1} ${y1} A ${r} ${r} 0 0 1 ${x2} ${y2}`;
                return (
                  <path
                    key={`excluded-arc-${idx}`}
                    d={d}
                    stroke="#ef4444"
                    strokeWidth="0.8"
                    fill="none"
                    opacity={0.7}
                  />
                );
              })}
            </g>
          )}
        </svg>

        {/* Four provider circles that push outward radially then fade */}
        {!hideIcons && (
        <AnimatePresence mode="popLayout">
          {activeItems.map((item) => {
            // Local per-item amplitude with delay (for fading only)
            const amp = localPhaseFor(nowMs, item.delayMs);
            // Outward then inward (to 50%) phase for radial motion
            const pushAmp = outwardInwardPhaseFor(nowMs, item.delayMs);
            // Use peak tick tip radius to prevent inward drift
            const tickPeakTipVB = outerRadius + LENGTH_BOOST + TICK_TRANSLATE_UNITS;
            // Fixed buffer (converted to viewBox units)
            const bufferVB = item.bufferPx * pxToVb;
            // Optional extra push to keep a touch of independent motion (no return)
            const extraPushVB = (item.pushPx * pxToVb) * pushAmp;
            // Optional static radial nudge per asset (px -> VB)
            const radialNudgeVB = (item.radialOffsetPx || 0) * pxToVb;
            // Constant base outward offset (px -> VB)
            const baseOffsetVB = BASE_CENTER_OFFSET_PX * pxToVb;
            // Center on the polar point with radial translation: peak tick tip + base offset + buffer + extra push + radial nudge
            const rCurrent = tickPeakTipVB + baseOffsetVB + bufferVB + extraPushVB + radialNudgeVB;
            const pos = posFromAngle(item.angleDeg, rCurrent);
            const dx = item.offsetPx?.dx ?? 0;
            const dy = item.offsetPx?.dy ?? 0;
            // Remove CSS translate centering; center by offsetting left/top by half size in VB units
            const halfChipVB = ((item.size + 15) / 2) * pxToVb;
            const left = `calc(${pos.left} - ${halfChipVB}%)`;
            const top = `calc(${pos.top} - ${halfChipVB}%)`;
            // Optional optical nudge in px -> VB
            const nudgeXVB = (dx || 0) * pxToVb;
            const nudgeYVB = (dy || 0) * pxToVb;
            const leftWithNudge = `calc(${left} + ${nudgeXVB}%)`;
            const topWithNudge = `calc(${top} + ${nudgeYVB}%)`;
            
          // Fade logic
          // - Before relax: fade in using local amplitude threshold
          // - During relax: fade out faster, fully gone by halfway through relax
          const relaxP = relaxProgressFor(nowMs, item.delayMs);
          const itemOpacity = relaxP > 0
            ? Math.max(0, 1 - Math.min(1, Math.pow(relaxP / 0.5, 1.2)))
            : (amp > FADE_OUT_START_AMP ? 1 : (amp / FADE_OUT_START_AMP));

          return (
          <motion.div
                key={item.id}
                className="absolute flex items-center justify-center"
            style={{ 
                  left: leftWithNudge,
                  top: topWithNudge,
                  pointerEvents: 'none'
            }}
                initial={{ opacity: 0, scale: 0.95 }}
                animate={{ opacity: itemOpacity, scale: 1 }}
                exit={{ opacity: 0, scale: 0.95 }}
            >
              <div
                className="rounded-full bg-white border border-gray-200 shadow-sm overflow-hidden flex items-center justify-center"
                style={{ 
                  width: item.size + 15, 
                  height: item.size + 15,
                  pointerEvents: 'none'
                }}
              >
                <img
                  src={item.src}
                  alt={item.alt}
                  className="object-cover block"
                  style={{ 
                    width: '100%', 
                    height: '100%',
                    transform: 'scale(1.25)',
                    pointerEvents: 'none'
                  }}
                  onError={() => { try { console.warn('[NullState] image failed', item); } catch (_) {} }}
                />
              </div>
            </motion.div>
            );
          })}
        </AnimatePresence>
        )}
        {DEBUG && !hideIcons && activeItems.map((item) => {
          const pushAmp = outwardInwardPhaseFor(nowMs, item.delayMs);
          const tickPeakTipVB = outerRadius + LENGTH_BOOST + TICK_TRANSLATE_UNITS;
          const bufferVB = item.bufferPx * pxToVb;
          const extraPushVB = (item.pushPx * pxToVb) * pushAmp;
          const radialNudgeVB = (item.radialOffsetPx || 0) * pxToVb;
          const baseOffsetVB = BASE_CENTER_OFFSET_PX * pxToVb;
          const rCurrent = tickPeakTipVB + baseOffsetVB + bufferVB + extraPushVB + radialNudgeVB;
          const pos = posFromAngle(item.angleDeg, rCurrent);
          const dotHalfVB = (4 / 2) * pxToVb;
          const left = `calc(${pos.left} - ${dotHalfVB}%)`;
          const top = `calc(${pos.top} - ${dotHalfVB}%)`;
          return (
            <div
              key={`dot-${item.id}`}
              className="absolute rounded-full"
              style={{
                left,
                top,
                width: 4,
                height: 4,
                background: 'red',
                pointerEvents: 'none'
              }}
            />
          );
        })}
      </div>
    </div>
  );
};

// Memoize to prevent unnecessary re-renders when parent re-renders
export const FloatingNullState = React.memo(FloatingNullStateComponent);

export default FloatingNullState;