// 1D Transfer Matrix Method calculator translated from the uploaded MATLAB Live Script.
// Units:
//   frequency inputs are displayed in THz but converted to Hz internally.
//   layer thicknesses are in meters.

const C0 = 299792458;

// -----------------------------------------------------------------------------
// Complex arithmetic
// -----------------------------------------------------------------------------
function c(re, im = 0) { return { re, im }; }
function cadd(a, b) { return c(a.re + b.re, a.im + b.im); }
function csub(a, b) { return c(a.re - b.re, a.im - b.im); }
function cmul(a, b) { return c(a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re); }
function cdiv(a, b) {
  const den = b.re * b.re + b.im * b.im;
  return c((a.re * b.re + a.im * b.im) / den, (a.im * b.re - a.re * b.im) / den);
}
function cabs2(a) { return a.re * a.re + a.im * a.im; }
function cabs(a) { return Math.sqrt(cabs2(a)); }
function cexp(a) {
  const er = Math.exp(a.re);
  return c(er * Math.cos(a.im), er * Math.sin(a.im));
}
function csqrt(z) {
  const r = Math.hypot(z.re, z.im);
  const re = Math.sqrt((r + z.re) / 2);
  const im = Math.sign(z.im || 1) * Math.sqrt(Math.max((r - z.re) / 2, 0));
  return c(re, im);
}

// -----------------------------------------------------------------------------
// 2x2 complex matrices
// -----------------------------------------------------------------------------
function m2(a00, a01, a10, a11) { return [[a00, a01], [a10, a11]]; }
function mmul(A, B) {
  return [
    [cadd(cmul(A[0][0], B[0][0]), cmul(A[0][1], B[1][0])),
     cadd(cmul(A[0][0], B[0][1]), cmul(A[0][1], B[1][1]))],
    [cadd(cmul(A[1][0], B[0][0]), cmul(A[1][1], B[1][0])),
     cadd(cmul(A[1][0], B[0][1]), cmul(A[1][1], B[1][1]))]
  ];
}
function mvecmul(A, v) {
  return [
    cadd(cmul(A[0][0], v[0]), cmul(A[0][1], v[1])),
    cadd(cmul(A[1][0], v[0]), cmul(A[1][1], v[1]))
  ];
}
function minv(A) {
  const det = csub(cmul(A[0][0], A[1][1]), cmul(A[0][1], A[1][0]));
  return [
    [cdiv(A[1][1], det), cdiv(cmul(c(-1), A[0][1]), det)],
    [cdiv(cmul(c(-1), A[1][0]), det), cdiv(A[0][0], det)]
  ];
}

// -----------------------------------------------------------------------------
// User-editable default structure from MATLAB script
// -----------------------------------------------------------------------------
const nAir = c(1, 0);
const nSi = c(3.42, 0);
const nDefect = c(3.42, 0);
const nGold = c(400, 500); // currently not used in default stack, preserved from MATLAB code.

const dAir = 197e-6;
const dSi = 50e-6;
const dDefect = 100e-6;
const dGold = 150e-9;

const nListDefault = [nSi, nAir, nSi, nAir, nDefect, nAir, nSi, nAir, nSi];
const dListDefault = [dSi, dAir, dSi, dAir, dDefect, dAir, dSi, dAir, dSi];

let lastResult = null;

// -----------------------------------------------------------------------------
// TMM core: equivalent to MATLAB TMM1D(n,d,freq)
// -----------------------------------------------------------------------------
function tmm1D(nList, dList, freqHz) {
  const k0 = 2 * Math.PI * freqHz / C0;
  const k = nList.map(n => cmul(n, c(k0, 0)));

  const D0 = m2(c(1), c(1), c(k0), c(-k0));
  const D = [];
  const Dinv = [];
  const P = [];
  const Pinv = [];

  let Mat = null;

  for (let l = 0; l < k.length; l++) {
    D[l] = m2(c(1), c(1), k[l], cmul(c(-1), k[l]));
    Dinv[l] = minv(D[l]);

    const phase = cmul(k[l], c(0, dList[l])); // i*k*d
    const expForward = cexp(phase);
    const expBackward = cexp(cmul(c(-1), phase));
    P[l] = m2(expForward, c(0), c(0), expBackward);
    Pinv[l] = minv(P[l]);

    const layerMat = mmul(mmul(D[l], Pinv[l]), Dinv[l]);

    if (l === 0) {
      Mat = mmul(minv(D0), layerMat);
    } else if (l === k.length - 1) {
      Mat = mmul(mmul(Mat, layerMat), D0);
    } else {
      Mat = mmul(Mat, layerMat);
    }
  }

  return { k, D, Dinv, P, Pinv, D0, Mat };
}

function calcTR(nList, dList, freqHz) {
  const { Mat } = tmm1D(nList, dList, freqHz);
  const T = 1 / cabs2(Mat[0][0]);
  const R = cabs2(cdiv(Mat[1][0], Mat[0][0]));
  return { T, R };
}

// Equivalent to MATLAB TMM1DEfield(n,d,f,pointsPerLayer) but returns arrays for plotting.
function tmm1DEfield(nList, dList, freqHz, samplesPerLayer = 500) {
  const { k, D, Dinv, P, D0, Mat } = tmm1D(nList, dList, freqHz);

  const vec0 = mvecmul(Mat, [c(1), c(0)]);
  const vecs = [];
  vecs[0] = mvecmul(mmul(Dinv[0], D0), vec0);

  const totalD = [];
  const totalE = [];
  const materialRegions = [];

  let offset = 0;

  for (let l = 0; l < nList.length; l++) {
    if (l > 0) {
      vecs[l] = mvecmul(mmul(mmul(Dinv[l], D[l - 1]), P[l - 1]), vecs[l - 1]);
    }

    if (Math.abs(nList[l].re - 1) > 1e-12 || Math.abs(nList[l].im) > 1e-12) {
      materialRegions.push({ x0: offset, x1: offset + dList[l] });
    }

    for (let i = 0; i < samplesPerLayer; i++) {
      const z = dList[l] * i / (samplesPerLayer - 1);
      const forward = cmul(vecs[l][0], cexp(cmul(k[l], c(0, z))));
      const backward = cmul(vecs[l][1], cexp(cmul(k[l], c(0, -z))));
      const E = cabs(cadd(forward, backward));
      totalD.push(offset + z);
      totalE.push(E);
    }

    offset += dList[l];
  }

  return { totalD, totalE, materialRegions, height: Math.max(...totalE) };
}

function findPeaks(y, x) {
  const peaks = [];
  for (let i = 1; i < y.length - 1; i++) {
    if (y[i] > y[i - 1] && y[i] >= y[i + 1]) {
      peaks.push({ index: i, x: x[i], y: y[i] });
    }
  }
  return peaks;
}

function readNumber(id) {
  return Number(document.getElementById(id).value);
}

function splitList(text) {
  return text
    .replace(/[;]+/g, ",")
    .split(/[\s,]+/)
    .map(v => v.trim())
    .filter(Boolean);
}

function parseComplexNumber(token) {
  let s = token.trim().toLowerCase().replace(/j/g, "i").replace(/\s+/g, "");
  if (!s) throw new Error("Empty complex number.");

  if (!s.includes("i")) {
    const re = Number(s);
    if (!Number.isFinite(re)) throw new Error(`Invalid refractive index: ${token}`);
    return c(re, 0);
  }

  s = s.replace(/i/g, "");
  if (s === "" || s === "+") return c(0, 1);
  if (s === "-") return c(0, -1);

  // Find the final + or - sign that separates real and imaginary parts.
  let splitIndex = -1;
  for (let i = 1; i < s.length; i++) {
    if ((s[i] === "+" || s[i] === "-") && s[i - 1] !== "e" && s[i - 1] !== "E") {
      splitIndex = i;
    }
  }

  if (splitIndex === -1) {
    const im = Number(s);
    if (!Number.isFinite(im)) throw new Error(`Invalid complex refractive index: ${token}`);
    return c(0, im);
  }

  const re = Number(s.slice(0, splitIndex));
  let imText = s.slice(splitIndex);
  if (imText === "+") imText = "1";
  if (imText === "-") imText = "-1";
  const im = Number(imText);

  if (!Number.isFinite(re) || !Number.isFinite(im)) {
    throw new Error(`Invalid complex refractive index: ${token}`);
  }
  return c(re, im);
}

function readLayerStructure() {
  const nTokens = splitList(document.getElementById("nListInput").value);
  const dTokens = splitList(document.getElementById("dListInput").value);

  if (nTokens.length === 0 || dTokens.length === 0) {
    throw new Error("Please enter both the refractive-index list and the thickness list.");
  }
  if (nTokens.length !== dTokens.length) {
    throw new Error(`The refractive-index list has ${nTokens.length} entries, but the thickness list has ${dTokens.length} entries. They must have the same length.`);
  }

  const nList = nTokens.map(parseComplexNumber);
  const dList = dTokens.map((token) => {
    const valueUm = Number(token);
    if (!Number.isFinite(valueUm) || valueUm <= 0) {
      throw new Error(`Invalid thickness: ${token}. Thickness values must be positive numbers in µm.`);
    }
    return valueUm * 1e-6;
  });

  return { nList, dList };
}

function linspaceFrequencyTHz(start, stop, step) {
  const f = [];
  const n = Math.floor((stop - start) / step + 1e-12) + 1;
  for (let i = 0; i < n; i++) {
    const val = start + i * step;
    if (val <= stop + step * 1e-9) f.push(val);
  }
  return f;
}

function runTMM() {
  const fStart = readNumber("fStart");
  const fStop = readNumber("fStop");
  const fStep = readNumber("fStep");

  if (!(fStart > 0) || !(fStop > fStart) || !(fStep > 0)) {
    alert("Please check that start frequency > 0, stop frequency > start frequency, and frequency step > 0.");
    return;
  }

  let layerStructure;
  try {
    layerStructure = readLayerStructure();
  } catch (err) {
    alert(err.message);
    return;
  }

  const { nList, dList } = layerStructure;
  const fTHz = linspaceFrequencyTHz(fStart, fStop, fStep);
  const T = [];
  const R = [];

  for (const f of fTHz) {
    const out = calcTR(nList, dList, f * 1e12);
    T.push(out.T);
    R.push(out.R);
  }

  lastResult = { fTHz, T, R, nList, dList };

  plotSpectra();
  updatePeaksAndField();
}

function plotSpectra() {
  if (!lastResult) return;

  const traceT = {
    x: lastResult.fTHz,
    y: lastResult.T,
    mode: "lines",
    name: "Transmittance, T"
  };

  const traceR = {
    x: lastResult.fTHz,
    y: lastResult.R,
    mode: "lines",
    name: "Reflectance, R"
  };

  Plotly.react("spectraPlot", [traceT, traceR], {
    margin: { l: 70, r: 25, t: 20, b: 65 },
    xaxis: { title: "Frequency (THz)" },
    yaxis: { title: "T / R" },
    legend: { orientation: "h", x: 0, y: 1.12 },
    hovermode: "x unified"
  }, { responsive: true });
}

function updatePeaksAndField() {
  if (!lastResult) {
    document.getElementById("peakList").textContent = "Run TMM to detect peaks.";
    document.getElementById("fieldInfo").textContent = "Run TMM to calculate the electric-field profile.";
    return;
  }

  const peakMin = readNumber("peakMin");
  const peakMax = readNumber("peakMax");
  const requestedPeakNumber = Math.max(1, Math.round(readNumber("peakNumber")));

  if (!(peakMax > peakMin)) {
    document.getElementById("peakList").textContent = "Freq max for peaks must be larger than Freq min for peaks.";
    document.getElementById("fieldInfo").textContent = "Electric-field profile not updated because the peak range is invalid.";
    Plotly.purge("fieldPlot");
    return;
  }

  const maxT = Math.max(...lastResult.T);
  const allPeaks = findPeaks(lastResult.T, lastResult.fTHz);
  const peaks = allPeaks.filter(p => p.y > 0.001 * maxT && p.x > peakMin && p.x < peakMax);

  const peakListText = peaks.length === 0
    ? "No peaks found in the selected frequency range."
    : peaks.map((p, i) => `${i + 1}: ${p.x.toFixed(6)} THz, T = ${p.y.toExponential(4)}`).join("\n");

  document.getElementById("peakList").textContent = peakListText;

  if (peaks.length === 0) {
    document.getElementById("fieldInfo").textContent = "No electric-field profile is shown because no peak was found in the selected range.";
    Plotly.purge("fieldPlot");
    return;
  }

  let peakNumber = requestedPeakNumber;
  if (peakNumber > peaks.length) {
    peakNumber = peaks.length;
    document.getElementById("peakNumber").value = peakNumber;
  }

  const selectedPeak = peaks[peakNumber - 1];
  plotElectricField(selectedPeak, peakNumber, peaks.length);
}

function plotElectricField(selectedPeak, peakNumber, totalPeaks) {
  const freqHz = selectedPeak.x * 1e12;
  const field = tmm1DEfield(lastResult.nList, lastResult.dList, freqHz, 500);

  const xUm = field.totalD.map(v => v * 1e6);

  const fieldTrace = {
    x: xUm,
    y: field.totalE,
    mode: "lines",
    name: "|E|"
  };

  const shapes = field.materialRegions.map(region => ({
    type: "rect",
    xref: "x",
    yref: "paper",
    x0: region.x0 * 1e6,
    x1: region.x1 * 1e6,
    y0: 0,
    y1: 1,
    fillcolor: "rgba(160,160,160,0.20)",
    line: { width: 0 },
    layer: "below"
  }));

  document.getElementById("fieldInfo").textContent =
    `Showing peak ${peakNumber} of ${totalPeaks}: frequency = ${selectedPeak.x.toFixed(6)} THz, T = ${selectedPeak.y.toExponential(4)}.`;

  Plotly.react("fieldPlot", [fieldTrace], {
    margin: { l: 70, r: 25, t: 20, b: 65 },
    xaxis: { title: "Position, d (µm)" },
    yaxis: { title: "Electric field, |E| (a.u.)" },
    shapes,
    hovermode: "x"
  }, { responsive: true });
}

// Event listeners
window.addEventListener("DOMContentLoaded", () => {
  document.getElementById("runButton").addEventListener("click", runTMM);

  // These controls do not need to recalculate the whole spectrum.
  // They only refilter the already-calculated peaks and then update the field profile.
  ["peakMin", "peakMax", "peakNumber"].forEach(id => {
    const el = document.getElementById(id);
    el.addEventListener("input", updatePeaksAndField);
    el.addEventListener("change", updatePeaksAndField);
  });

  runTMM();
});
