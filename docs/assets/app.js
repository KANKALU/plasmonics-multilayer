/* Reflectance explorer.
 *
 * Curves come from docs/data/curves.json, written by
 * scripts/export_web_data.py, which runs the same transfer-matrix code
 * as the rest of the repository. Nothing is recomputed here; the page
 * only draws.
 */

const SVG_NS = "http://www.w3.org/2000/svg";

const LAYER_FILL = {
  gold: "var(--au)",
  copper: "var(--cu)",
  ta2o5: "var(--ta)",
  sio2: "var(--ox)",
  quartz: "var(--sub)",
};

const MATERIAL_NAME = {
  gold: "Au", copper: "Cu", ta2o5: "Ta₂O₅", sio2: "SiO₂", quartz: "SiO₂ sub.", bk7: "BK7", air: "air",
};

/* Critical angle of the BK7/air interface: everything below it is
   ordinary refraction, everything above is total internal reflection. */
const THETA_C = 41.0;
/* Above this the immersion oil dominates the measurement. */
const OIL_ONSET = 70.0;

const state = { data: null, sample: "M7", lambda: "700", theta: 55, lang: "es" };

const T = {
  noMinimum: { es: "Sin mínimo en TE.", en: "No TE minimum." },
  teMin: {
    es: (r, t) => `Mínimo en TE: R = ${r.toFixed(3)} a ${t.toFixed(1)}°.`,
    en: (r, t) => `TE minimum: R = ${r.toFixed(3)} at ${t.toFixed(1)}°.`,
  },
  tmMin: {
    es: (r, t) => `Mínimo en TM: R = ${r.toFixed(3)} a ${t.toFixed(1)}°.`,
    en: (r, t) => `TM minimum: R = ${r.toFixed(3)} at ${t.toFixed(1)}°.`,
  },
  teYes: { es: "mínimo TE", en: "TE minimum" },
  teNo: { es: "sin mínimo TE", en: "no TE minimum" },
  angleAxis: { es: "Ángulo de incidencia θ (°)", en: "Angle of incidence θ (°)" },
  reflAxis: { es: "Reflectancia", en: "Reflectance" },
  oilBand: { es: "aceite", en: "oil" },
  loadFailed: {
    es: "No se pudieron cargar las curvas. Sirve el sitio por HTTP, no abriendo el archivo directamente.",
    en: "Could not load the curves. Serve the site over HTTP rather than opening the file directly.",
  },
};

const el = (tag, attrs = {}, parent = null) => {
  const node = document.createElementNS(SVG_NS, tag);
  for (const [k, v] of Object.entries(attrs)) node.setAttribute(k, v);
  if (parent) parent.appendChild(node);
  return node;
};

const clear = (node) => { while (node.firstChild) node.removeChild(node.firstChild); };

/* ------------------------------------------------------------------ */
/* language                                                            */
/* ------------------------------------------------------------------ */

function applyLanguage(lang) {
  state.lang = lang;
  document.documentElement.lang = lang;
  document.querySelectorAll("[data-es][data-en]").forEach((node) => {
    node.textContent = node.dataset[lang];
  });
  document.querySelectorAll(".lang-btn").forEach((btn) => {
    btn.classList.toggle("is-on", btn.dataset.lang === lang);
  });
  if (state.data) { drawChart(); drawSampleGrid(); updateReadout(); }
}

/* ------------------------------------------------------------------ */
/* geometry panel                                                      */
/* ------------------------------------------------------------------ */

function drawGeometry() {
  const svg = document.getElementById("geometry");
  clear(svg);
  const sample = state.data.samples[state.sample];

  const W = 440, H = 320;
  const surfaceY = 190;          // the prism base / stack top boundary
  const cx = 190;                // where the beam meets the surface

  /* prism: a hemicylindrical lens sitting on the sample, flat side down */
  const prismR = 150;
  el("path", {
    d: `M ${cx - prismR} ${surfaceY} A ${prismR} ${prismR} 0 0 1 ${cx + prismR} ${surfaceY} Z`,
    fill: "rgba(255,255,255,.72)", stroke: "var(--rule)", "stroke-width": 1.2,
  }, svg);
  el("text", {
    x: 62, y: surfaceY - 12,
    "font-family": "var(--mono)", "font-size": 11, fill: "var(--ink-3)",
  }, svg).textContent = "BK7";

  /* the deposited stack, drawn downward, to scale among themselves.
     The 500 nm substrate is compressed so 20 nm films stay visible. */
  const media = sample.media.slice(1, -1);
  const nm = sample.layers_nm;
  const drawn = media.map((m, i) => (m === "quartz" ? nm[i] * 0.045 : nm[i]));
  const total = drawn.reduce((a, b) => a + b, 0);
  const stackPx = 96;
  const scale = stackPx / total;

  let y = surfaceY;
  for (let i = media.length - 1; i >= 0; i -= 1) {
    const h = drawn[i] * scale;
    el("rect", {
      x: 40, y, width: 300, height: h,
      fill: LAYER_FILL[media[i]] || "var(--rule)",
      stroke: "rgba(20,24,28,.28)", "stroke-width": .6,
    }, svg);
    if (h > 9) {
      el("text", {
        x: 348, y: y + h / 2 + 3.5,
        "font-family": "var(--mono)", "font-size": 10, fill: "var(--ink-2)",
      }, svg).textContent = `${MATERIAL_NAME[media[i]] || media[i]} ${nm[i]}`;
    }
    y += h;
  }

  /* beams: incident and reflected, mirrored about the surface normal */
  const rad = (state.theta * Math.PI) / 180;
  const L = 132;
  const dx = L * Math.sin(rad), dy = L * Math.cos(rad);

  el("line", {
    x1: cx, y1: surfaceY, x2: cx, y2: surfaceY - 118,
    stroke: "var(--ink-3)", "stroke-width": .8, "stroke-dasharray": "3 3",
  }, svg);

  const beam = (x2, color) => el("line", {
    x1: cx, y1: surfaceY, x2, y2: surfaceY - dy,
    stroke: color, "stroke-width": 2.1, "stroke-linecap": "round",
  }, svg);
  beam(cx - dx, "var(--ink)");
  beam(cx + dx, "var(--te)");

  /* angle arc between the normal and the incident beam */
  const r = 34;
  const ax = cx - r * Math.sin(rad), ay = surfaceY - r * Math.cos(rad);
  el("path", {
    d: `M ${cx} ${surfaceY - r} A ${r} ${r} 0 0 0 ${ax} ${ay}`,
    fill: "none", stroke: "var(--ink-3)", "stroke-width": 1,
  }, svg);
  el("text", {
    x: cx - (r + 13) * Math.sin(rad / 2), y: surfaceY - (r + 13) * Math.cos(rad / 2) + 4,
    "text-anchor": "middle", "font-family": "var(--mono)", "font-size": 11, fill: "var(--ink-2)",
  }, svg).textContent = "θ";

  el("text", {
    x: 40, y: H - 14, "font-family": "var(--mono)", "font-size": 10, fill: "var(--ink-3)",
  }, svg).textContent = `${state.sample} · ${sample.label}`;
}

/* ------------------------------------------------------------------ */
/* chart panel                                                         */
/* ------------------------------------------------------------------ */

function drawChart() {
  const svg = document.getElementById("chart");
  clear(svg);

  const W = 520, H = 320;
  const m = { top: 18, right: 16, bottom: 44, left: 52 };
  const iw = W - m.left - m.right;
  const ih = H - m.top - m.bottom;

  const theta = state.data.theta_deg;
  const curves = state.data.samples[state.sample].curves[state.lambda];
  const t0 = theta[0], t1 = theta[theta.length - 1];

  const X = (t) => m.left + ((t - t0) / (t1 - t0)) * iw;
  const Y = (r) => m.top + (1 - r) * ih;

  /* the region the report does not interpret, shaded rather than hidden */
  el("rect", {
    x: X(OIL_ONSET), y: m.top, width: X(t1) - X(OIL_ONSET), height: ih,
    fill: "rgba(20,24,28,.05)",
  }, svg);
  el("text", {
    x: X(OIL_ONSET) + 5, y: m.top + 12,
    "font-family": "var(--mono)", "font-size": 9, fill: "var(--ink-3)",
  }, svg).textContent = T.oilBand[state.lang];

  /* grid and axes */
  for (let r = 0; r <= 1.0001; r += 0.25) {
    el("line", {
      x1: m.left, y1: Y(r), x2: m.left + iw, y2: Y(r),
      stroke: "var(--rule-soft)", "stroke-width": 1,
    }, svg);
    el("text", {
      x: m.left - 8, y: Y(r) + 3.5, "text-anchor": "end",
      "font-family": "var(--mono)", "font-size": 10, fill: "var(--ink-3)",
    }, svg).textContent = r.toFixed(2);
  }
  for (let t = 20; t <= 80; t += 10) {
    el("line", {
      x1: X(t), y1: m.top + ih, x2: X(t), y2: m.top + ih + 4,
      stroke: "var(--rule)", "stroke-width": 1,
    }, svg);
    el("text", {
      x: X(t), y: m.top + ih + 17, "text-anchor": "middle",
      "font-family": "var(--mono)", "font-size": 10, fill: "var(--ink-3)",
    }, svg).textContent = t;
  }

  /* critical angle */
  el("line", {
    x1: X(THETA_C), y1: m.top, x2: X(THETA_C), y2: m.top + ih,
    stroke: "var(--ink-3)", "stroke-width": 1, "stroke-dasharray": "4 3",
  }, svg);
  el("text", {
    x: X(THETA_C) + 5, y: m.top + ih - 6,
    "font-family": "var(--mono)", "font-size": 9, fill: "var(--ink-3)",
  }, svg).textContent = "θc";

  const path = (values, color, width) => {
    const d = values.map((r, i) => `${i ? "L" : "M"} ${X(theta[i]).toFixed(2)} ${Y(r).toFixed(2)}`).join(" ");
    el("path", { d, fill: "none", stroke: color, "stroke-width": width, "stroke-linejoin": "round" }, svg);
  };
  path(curves.p, "var(--tm)", 2);
  path(curves.s, "var(--te)", 2.4);

  /* current angle marker */
  const i = nearestIndex(theta, state.theta);
  el("line", {
    x1: X(theta[i]), y1: m.top, x2: X(theta[i]), y2: m.top + ih,
    stroke: "var(--ink)", "stroke-width": 1, opacity: .38,
  }, svg);
  el("circle", { cx: X(theta[i]), cy: Y(curves.p[i]), r: 4, fill: "var(--tm)", stroke: "var(--panel)", "stroke-width": 1.5 }, svg);
  el("circle", { cx: X(theta[i]), cy: Y(curves.s[i]), r: 4, fill: "var(--te)", stroke: "var(--panel)", "stroke-width": 1.5 }, svg);

  /* axis titles */
  el("text", {
    x: m.left + iw / 2, y: H - 8, "text-anchor": "middle",
    "font-family": "var(--sans)", "font-size": 11, fill: "var(--ink-2)",
  }, svg).textContent = T.angleAxis[state.lang];
  el("text", {
    x: 13, y: m.top + ih / 2, "text-anchor": "middle",
    transform: `rotate(-90 13 ${m.top + ih / 2})`,
    "font-family": "var(--sans)", "font-size": 11, fill: "var(--ink-2)",
  }, svg).textContent = T.reflAxis[state.lang];

  /* inline legend, so the curves need no external key */
  const legend = [["TM (p)", "var(--tm)"], ["TE (s)", "var(--te)"]];
  legend.forEach(([label, color], k) => {
    const lx = m.left + 10, ly = m.top + 14 + k * 16;
    el("line", { x1: lx, y1: ly, x2: lx + 16, y2: ly, stroke: color, "stroke-width": 2.4 }, svg);
    el("text", {
      x: lx + 22, y: ly + 3.5,
      "font-family": "var(--mono)", "font-size": 10, fill: "var(--ink-2)",
    }, svg).textContent = label;
  });
}

/* ------------------------------------------------------------------ */
/* readout                                                             */
/* ------------------------------------------------------------------ */

const nearestIndex = (arr, value) => {
  let best = 0, bestD = Infinity;
  for (let i = 0; i < arr.length; i += 1) {
    const d = Math.abs(arr[i] - value);
    if (d < bestD) { bestD = d; best = i; }
  }
  return best;
};

/* Deepest point strictly inside the interpretable window. A result at
   the very edge means the curve is monotonic there, i.e. no minimum. */
function interiorMinimum(theta, values, lo = THETA_C + 1, hi = OIL_ONSET) {
  let idx = -1, best = Infinity;
  for (let i = 0; i < theta.length; i += 1) {
    if (theta[i] < lo || theta[i] > hi) continue;
    if (values[i] < best) { best = values[i]; idx = i; }
  }
  if (idx <= 0 || idx >= theta.length - 1) return null;
  const isTurning = values[idx] < values[idx - 1] && values[idx] < values[idx + 1];
  if (!isTurning || theta[idx] <= lo + 0.5 || theta[idx] >= hi - 0.5) return null;
  return { r: values[idx], theta: theta[idx] };
}

function updateReadout() {
  const theta = state.data.theta_deg;
  const curves = state.data.samples[state.sample].curves[state.lambda];
  const i = nearestIndex(theta, state.theta);

  document.getElementById("r-tm").textContent = curves.p[i].toFixed(3);
  document.getElementById("r-te").textContent = curves.s[i].toFixed(3);
  document.getElementById("theta-out").textContent = `${state.theta.toFixed(2)}°`;

  const te = interiorMinimum(theta, curves.s);
  const tm = interiorMinimum(theta, curves.p);
  const parts = [
    te ? T.teMin[state.lang](te.r, te.theta) : T.noMinimum[state.lang],
    tm ? T.tmMin[state.lang](tm.r, tm.theta) : null,
  ].filter(Boolean);
  document.getElementById("sample-note").textContent = parts.join(" ");
}

/* ------------------------------------------------------------------ */
/* sample cards                                                        */
/* ------------------------------------------------------------------ */

function stackThumb(sample) {
  const svg = document.createElementNS(SVG_NS, "svg");
  svg.setAttribute("viewBox", "0 0 160 92");
  svg.setAttribute("role", "img");

  const media = sample.media.slice(1, -1);
  const nm = sample.layers_nm;
  const drawn = media.map((m, i) => (m === "quartz" ? nm[i] * 0.045 : nm[i]));
  const total = drawn.reduce((a, b) => a + b, 0);
  const scale = 78 / total;

  let y = 92 - 78;
  for (let i = media.length - 1; i >= 0; i -= 1) {
    const h = drawn[i] * scale;
    el("rect", {
      x: 0, y, width: 160, height: h,
      fill: LAYER_FILL[media[i]] || "var(--rule)",
      stroke: "rgba(20,24,28,.25)", "stroke-width": .5,
    }, svg);
    y += h;
  }
  return svg;
}

function drawSampleGrid() {
  const grid = document.getElementById("sample-grid");
  grid.innerHTML = "";

  for (const [key, sample] of Object.entries(state.data.samples)) {
    const hasTE = Object.values(sample.curves).some(
      (c) => interiorMinimum(state.data.theta_deg, c.s) !== null
    );

    const card = document.createElement("article");
    card.className = "sample-card";

    const h = document.createElement("h4");
    h.textContent = key;
    card.appendChild(h);
    card.appendChild(stackThumb(sample));

    const label = document.createElement("p");
    label.className = "stack-label";
    label.textContent = sample.label;
    card.appendChild(label);

    const flag = document.createElement("span");
    flag.className = `te-flag ${hasTE ? "te-yes" : "te-no"}`;
    flag.textContent = hasTE ? T.teYes[state.lang] : T.teNo[state.lang];
    card.appendChild(flag);

    grid.appendChild(card);
  }
}

/* ------------------------------------------------------------------ */
/* controls and boot                                                   */
/* ------------------------------------------------------------------ */

function buildChips(containerId, values, current, onPick) {
  const box = document.getElementById(containerId);
  box.innerHTML = "";
  values.forEach((value) => {
    const btn = document.createElement("button");
    btn.type = "button";
    btn.className = `chip${String(value) === String(current) ? " is-on" : ""}`;
    btn.textContent = value;
    btn.addEventListener("click", () => {
      box.querySelectorAll(".chip").forEach((c) => c.classList.remove("is-on"));
      btn.classList.add("is-on");
      onPick(String(value));
    });
    box.appendChild(btn);
  });
}

function redraw() { drawGeometry(); drawChart(); updateReadout(); }

async function boot() {
  document.querySelectorAll(".lang-btn").forEach((btn) => {
    btn.addEventListener("click", () => applyLanguage(btn.dataset.lang));
  });
  applyLanguage("es");

  let data;
  try {
    const res = await fetch("data/curves.json");
    if (!res.ok) throw new Error(res.status);
    data = await res.json();
  } catch (err) {
    document.getElementById("sample-note").textContent = T.loadFailed[state.lang];
    return;
  }
  state.data = data;

  buildChips("sample-chips", Object.keys(data.samples), state.sample, (v) => {
    state.sample = v;
    redraw();
  });
  buildChips("lambda-chips", data.wavelengths_nm.map((w) => `${w}`), state.lambda, (v) => {
    state.lambda = v;
    redraw();
  });

  const slider = document.getElementById("theta");
  slider.addEventListener("input", () => {
    state.theta = parseFloat(slider.value);
    redraw();
  });

  drawSampleGrid();
  redraw();
}

boot();
