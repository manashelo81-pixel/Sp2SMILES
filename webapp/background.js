/* Interactive page background for the Streamlit app.
 * MODE: "electrons" (electrons spiralling around a nucleus that follows the mouse)
 *       "spectrum"  (ions arriving and building mass-spectrum peaks; move the mouse to scan m/z, click to add ions)
 *       "off"
 * The script runs inside a Streamlit component iframe and draws on a canvas it adds to the parent page. */
(function (MODE) {
  var P = window.parent, D = P.document;
  if (P.__bgStop) { try { P.__bgStop(); } catch (e) {} P.__bgStop = null; }
  var old = D.getElementById('bg-canvas'); if (old) old.remove();
  var styleEl = D.getElementById('bg-style');
  if (MODE === 'off') { if (styleEl) styleEl.remove(); return; }

  /* ---------- theme (follows the Streamlit light/dark setting) ---------- */
  var app = D.querySelector('.stApp');
  var rgb = ((P.getComputedStyle(app || D.body).color || '').match(/[\d.]+/g) || [30, 30, 30]).map(Number);
  var dark = (0.299 * rgb[0] + 0.587 * rgb[1] + 0.114 * rgb[2]) > 150;
  var css =
    'html,body{min-height:100vh}' +
    'body{background-image:' + (dark ? 'linear-gradient(135deg,#0a1122,#101b38 60%,#0c1529)'
                                      : 'linear-gradient(135deg,#f7faff,#e5edfb 60%,#f4f8ff)') +
      ' !important;background-attachment:fixed !important;background-size:cover !important;background-repeat:no-repeat !important}' +
    '.stApp,[data-testid="stAppViewContainer"],[data-testid="stMain"],[data-testid="stHeader"],[data-testid="stToolbar"]{background:transparent !important}' +
        '[data-testid="stSidebar"]{background:' + (dark ? 'rgba(12,17,32,.55)' : 'rgba(255,255,255,.55)') + ' !important}' +
    '[data-testid="stMainBlockContainer"],.block-container{background:' + (dark ? 'rgba(12,17,32,.42)' : 'rgba(255,255,255,.5)') + ';border-radius:16px}';
  if (!styleEl) { styleEl = D.createElement('style'); styleEl.id = 'bg-style'; D.head.appendChild(styleEl); }
  styleEl.textContent = css;

  /* ---------- canvas ---------- */
  var cv = D.createElement('canvas');
  cv.id = 'bg-canvas';
  cv.style.cssText = 'position:fixed;top:0;left:0;width:100vw;height:100vh;z-index:-1;pointer-events:none';
  D.body.appendChild(cv);
  var ctx = cv.getContext('2d');
  var W = 0, H = 0;
  function resize() {
    var dpr = Math.min(P.devicePixelRatio || 1, 2);
    W = P.innerWidth; H = P.innerHeight;
    cv.width = Math.round(W * dpr); cv.height = Math.round(H * dpr);
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  }
  resize();
  var COL = dark ? ['#7cc4ff', '#a78bfa', '#5eead4', '#fcd34d'] : ['#2563eb', '#7c3aed', '#0d9488', '#ea580c'];
  var INK = dark ? '230,240,255' : '30,50,90';
  function rnd(a, b) { return a + Math.random() * (b - a); }

  /* ---------- pointer ---------- */
  var mouse = { x: W / 2, y: H / 2, seen: -1e9, down: 0 };
  function onMove(e) {
    var t = e.touches && e.touches[0] ? e.touches[0] : e;
    mouse.x = t.clientX; mouse.y = t.clientY; mouse.seen = P.performance ? P.performance.now() : Date.now();
  }
  function onClick(e) { onMove(e); mouse.down = 1; }
  P.addEventListener('mousemove', onMove, { passive: true });
  P.addEventListener('touchmove', onMove, { passive: true });
  P.addEventListener('click', onClick, { passive: true });
  P.addEventListener('resize', resize);

  /* ======================= ELECTRONS ======================= */
  var E = [], cx = W / 2, cy = H / 2, TRAIL = 30, ARMS = 4;
  function spawn(edge) {
    var maxR = Math.hypot(W, H) * 0.5;
    var r = edge ? rnd(0.6, 1) * maxR : 25 + Math.pow(Math.random(), 1.6) * (maxR - 25);   // denser near the nucleus
    var arm = ((Math.random() * ARMS) | 0) * (6.283 / ARMS);
    return { a: arm + Math.log(r) * 1.7 + rnd(-0.22, 0.22), r: r,            // electrons start along spiral arms
             dir: Math.random() < 0.9 ? 1 : -1, k: rnd(0.75, 1.35), c: COL[(Math.random() * COL.length) | 0],
             s: rnd(1.8, 3.2), trail: [] };
  }
  function initElectrons() {
    E = [];
    var n = Math.round(Math.max(90, Math.min(190, W * H / 6500)));
    for (var i = 0; i < n; i++) E.push(spawn(false));
  }
  function stepElectrons(dt, t) {
    var idle = (t - mouse.seen) > 4000;
    var tx = idle ? W * (0.5 + 0.28 * Math.sin(t / 5200)) : mouse.x;
    var ty = idle ? H * (0.5 + 0.22 * Math.sin(t / 3700 + 1.3)) : mouse.y;
    cx += (tx - cx) * Math.min(1, dt * 3); cy += (ty - cy) * Math.min(1, dt * 3);
    if (mouse.down) {                                    // click: release a burst of electrons from the nucleus
      for (var b = 0; b < 14; b++) { var s = spawn(false); s.r = rnd(10, 60); E.push(s); }
      if (E.length > 300) E.splice(0, E.length - 300);
      mouse.down = 0;
    }
    ctx.clearRect(0, 0, W, H);
    var g = ctx.createRadialGradient(cx, cy, 0, cx, cy, 130);       // nucleus glow
    g.addColorStop(0, 'rgba(' + (dark ? '124,196,255' : '37,99,235') + ',0.38)');
    g.addColorStop(1, 'rgba(' + (dark ? '124,196,255' : '37,99,235') + ',0)');
    ctx.fillStyle = g; ctx.beginPath(); ctx.arc(cx, cy, 130, 0, 6.283); ctx.fill();
    ctx.fillStyle = dark ? 'rgba(180,220,255,.9)' : 'rgba(37,99,235,.9)';
    ctx.beginPath(); ctx.arc(cx, cy, 3.2, 0, 6.283); ctx.fill();
    ctx.lineCap = 'round';
    for (var i = 0; i < E.length; i++) {
      var p = E[i];
      p.a += p.dir * p.k * (150 / (p.r + 60)) * dt;                // inner electrons orbit faster
      p.r -= (14 + p.r * 0.07) * p.k * dt;                         // slow inward spiral
      if (p.r < 9) { E[i] = p = spawn(true); }
      var x = cx + p.r * Math.cos(p.a), y = cy + p.r * Math.sin(p.a);
      p.trail.push(x, y);
      if (p.trail.length > TRAIL * 2) p.trail.splice(0, 2);
      var tr = p.trail, m = tr.length / 2;
      ctx.strokeStyle = p.c; ctx.lineWidth = p.s * 1.0;
      for (var j = 1; j < m; j++) {
        ctx.globalAlpha = 0.6 * j / m;
        ctx.beginPath(); ctx.moveTo(tr[2 * j - 2], tr[2 * j - 1]); ctx.lineTo(tr[2 * j], tr[2 * j + 1]); ctx.stroke();
      }
      ctx.globalAlpha = 0.9; ctx.fillStyle = p.c;
      ctx.beginPath(); ctx.arc(x, y, p.s, 0, 6.283); ctx.fill();
    }
    ctx.globalAlpha = 1;
  }

  /* ======================= SPECTRUM ======================= */
  var PK = [], IONS = [], rate = 0, BASE = 34;
  function addPeak(x) {
    var w = Math.pow(Math.random(), 2.2) * 0.95 + 0.05;                  // a few dominant peaks, many small ones
    PK.push({ x: x, w: w, target: (0.05 + 0.27 * w) * H, h: 0, flash: 0, c: COL[(Math.random() * 2) | 0] });
    PK[PK.length - 1].h = PK[PK.length - 1].target * rnd(0.3, 0.8);
  }
  function pickPeak() {
    var tot = 0, i; for (i = 0; i < PK.length; i++) tot += PK[i].w;
    var r = Math.random() * tot; for (i = 0; i < PK.length; i++) { r -= PK[i].w; if (r <= 0) return PK[i]; }
    return PK[PK.length - 1];
  }
  function initSpectrum() {
    PK = []; IONS = [];
    var n = Math.round(Math.max(18, Math.min(52, W / 34)));
    for (var i = 0; i < n; i++) addPeak(rnd(0.03, 0.97) * W);
    PK.sort(function (a, b) { return a.x - b.x; });
  }
  function drop(x, pk) { IONS.push({ x: x + rnd(-3, 3), y: -8, vy: rnd(120, 260), pk: pk, c: pk ? pk.c : COL[0] }); }
  function stepSpectrum(dt, t) {
    var y0 = H - BASE, near = (t - mouse.seen) < 2500;
    ctx.clearRect(0, 0, W, H);
    /* ions arrive at random peaks, and quickly where the mouse points ("scanning" m/z) */
    rate += dt * 38;
    while (rate >= 1) { rate -= 1; if (PK.length) { var q = pickPeak(); drop(q.x, q); } }
    if (near && Math.random() < dt * 28) {
      var best = null, bd = 1e9;
      for (var k = 0; k < PK.length; k++) { var d = Math.abs(PK[k].x - mouse.x); if (d < bd) { bd = d; best = PK[k]; } }
      if (best && bd < 60) drop(best.x, best); else drop(mouse.x, null);
    }
    if (mouse.down) {                                       // click: new peak (or a burst on a nearby one)
      var hit = null;
      for (var k2 = 0; k2 < PK.length; k2++) if (Math.abs(PK[k2].x - mouse.x) < 14) hit = PK[k2];
      if (!hit && PK.length < 70) { addPeak(mouse.x); hit = PK[PK.length - 1]; }
      if (hit) for (var b = 0; b < 25; b++) { var ion = { x: hit.x + rnd(-2, 2), y: -rnd(0, 200), vy: rnd(160, 300), pk: hit, c: hit.c }; IONS.push(ion); }
      mouse.down = 0;
    }
    /* baseline + m/z axis */
    ctx.strokeStyle = 'rgba(' + INK + ',0.22)'; ctx.lineWidth = 1;
    ctx.beginPath(); ctx.moveTo(0, y0 + 0.5); ctx.lineTo(W, y0 + 0.5); ctx.stroke();
    ctx.fillStyle = 'rgba(' + INK + ',0.35)'; ctx.font = '10px sans-serif'; ctx.textAlign = 'center';
    for (var mz = 100; mz < 1000; mz += 100) {
      var tx = mz / 1000 * W; ctx.beginPath(); ctx.moveTo(tx + 0.5, y0); ctx.lineTo(tx + 0.5, y0 + 5); ctx.stroke();
      ctx.fillText(String(mz), tx, y0 + 17);
    }
    /* peaks */
    ctx.lineCap = 'round';
    for (var i = 0; i < PK.length; i++) {
      var p = PK[i];
      p.h *= (1 - 0.09 * dt);                               // slow decay; arriving ions keep rebuilding the peak
      p.flash *= (1 - 5 * dt);
      var top = y0 - p.h;
      ctx.strokeStyle = p.c; ctx.globalAlpha = 0.75 + 0.25 * Math.min(1, p.flash);
      ctx.lineWidth = 3; ctx.beginPath(); ctx.moveTo(p.x, y0); ctx.lineTo(p.x, top); ctx.stroke();
      if (p.h > 8) {
        ctx.fillStyle = p.c; ctx.beginPath(); ctx.arc(p.x, top, 2.6 + 2.5 * p.flash, 0, 6.283); ctx.fill();
        if (p.h > H * 0.2) { ctx.globalAlpha = 0.5; ctx.fillStyle = 'rgba(' + INK + ',1)'; ctx.font = '10px sans-serif';
          ctx.textAlign = 'center'; ctx.fillText((p.x / W * 1000).toFixed(1), p.x, top - 8); }
      }
    }
    /* falling ions */
    for (var n = IONS.length - 1; n >= 0; n--) {
      var io = IONS[n];
      io.vy += 340 * dt; io.y += io.vy * dt;
      var stop = io.pk ? y0 - io.pk.h : y0;
      if (io.y >= stop) {
        if (io.pk) { io.pk.h = Math.min(io.pk.h + io.pk.target * rnd(0.04, 0.09), io.pk.target * 1.15); io.pk.flash = 1; }
        IONS.splice(n, 1); continue;
      }
      ctx.globalAlpha = 0.85; ctx.fillStyle = io.c;
      ctx.beginPath(); ctx.arc(io.x, io.y, 2, 0, 6.283); ctx.fill();
      ctx.globalAlpha = 0.25; ctx.fillRect(io.x - 0.7, io.y - 14, 1.4, 14);
    }
    /* cursor read-out */
    if (near) {
      ctx.globalAlpha = 0.28; ctx.strokeStyle = 'rgba(' + INK + ',1)'; ctx.setLineDash([3, 4]);
      ctx.beginPath(); ctx.moveTo(mouse.x + 0.5, y0 - H * 0.45); ctx.lineTo(mouse.x + 0.5, y0); ctx.stroke(); ctx.setLineDash([]);
      ctx.globalAlpha = 0.6; ctx.fillStyle = 'rgba(' + INK + ',1)'; ctx.font = '11px sans-serif'; ctx.textAlign = 'center';
      ctx.fillText('m/z ' + (mouse.x / W * 1000).toFixed(1), mouse.x, y0 - H * 0.45 - 6);
    }
    ctx.globalAlpha = 1;
  }

  /* ---------- safety net: undo all changes if the main content area looks collapsed ---------- */
  function revert() {
    if (P.__bgStop) { try { P.__bgStop(); } catch (e) {} P.__bgStop = null; }
    var s = D.getElementById('bg-style'); if (s) s.remove();
  }
  function check() {
    var mc = D.querySelector('[data-testid="stMainBlockContainer"]') || D.querySelector('.block-container');
    if (mc && mc.getBoundingClientRect) {
      var r = mc.getBoundingClientRect();
      if (r.height < 40 || r.width < 40) revert();
    }
  }
  setTimeout(check, 1200); setTimeout(check, 4000);

  /* ---------- loop ---------- */
  if (MODE === 'spectrum') initSpectrum(); else initElectrons();
  var raf = 0, last = 0, stopped = false;
  function frame(t) {
    if (stopped) return;
    var dt = last ? Math.min(0.05, (t - last) / 1000) : 0.016; last = t;
    if (MODE === 'spectrum') stepSpectrum(dt, t); else stepElectrons(dt, t);
    raf = P.requestAnimationFrame(frame);
  }
  raf = P.requestAnimationFrame(frame);
  P.__bgStop = function () {
    stopped = true; P.cancelAnimationFrame(raf);
    P.removeEventListener('mousemove', onMove); P.removeEventListener('touchmove', onMove);
    P.removeEventListener('click', onClick); P.removeEventListener('resize', resize);
    var c = D.getElementById('bg-canvas'); if (c) c.remove();
  };
})(__MODE__);
