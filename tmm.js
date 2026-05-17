/*
  1D Transfer Matrix Method calculator converted from MATLAB Live Script.
  Original MATLAB logic:
    - TMM1D(n, d, freq)
    - TMM1DEfield(n, d, f)
    - optional quantum-well cyclotron-resonance dielectric function

  Units:
    - frequency: Hz
    - thickness: m
    - n: complex refractive index
*/

const C0 = 299792458;
const EPS0 = 8.85418782e-12;
const QE = 1.60217662e-19;
const ME = 9.11e-31;

class Complex {
  constructor(re = 0, im = 0) {
    this.re = re;
    this.im = im;
  }
  static from(x) {
    if (x instanceof Complex) return x;
    if (typeof x === "number") return new Complex(x, 0);
    if (typeof x === "object" && x !== null && "re" in x) return new Complex(x.re, x.im || 0);
    throw new Error("Cannot convert to Complex: " + x);
  }
  add(z) { z = Complex.from(z); return new Complex(this.re + z.re, this.im + z.im); }
  sub(z) { z = Complex.from(z); return new Complex(this.re - z.re, this.im - z.im); }
  mul(z) {
    z = Complex.from(z);
    return new Complex(this.re * z.re - this.im * z.im, this.re * z.im + this.im * z.re);
  }
  div(z) {
    z = Complex.from(z);
    const den = z.re * z.re + z.im * z.im;
    return new Complex((this.re * z.re + this.im * z.im) / den, (this.im * z.re - this.re * z.im) / den);
  }
  neg() { return new Complex(-this.re, -this.im); }
  abs() { return Math.hypot(this.re, this.im); }
  exp() {
    const e = Math.exp(this.re);
    return new Complex(e * Math.cos(this.im), e * Math.sin(this.im));
  }
  sqrt() {
    const r = this.abs();
    const theta = Math.atan2(this.im, this.re) / 2;
    return new Complex(Math.sqrt(r) * Math.cos(theta), Math.sqrt(r) * Math.sin(theta));
  }
  toString(ndigits = 4) {
    const a = this.re.toFixed(ndigits);
    const b = Math.abs(this.im).toFixed(ndigits);
    const s = this.im >= 0 ? "+" : "-";
    return `${a} ${s} ${b}i`;
  }
}

const I = new Complex(0, 1);
const C = (re, im = 0) => new Complex(re, im);

function cadd(a, b) { return Complex.from(a).add(b); }
function csub(a, b) { return Complex.from(a).sub(b); }
function cmul(a, b) { return Complex.from(a).mul(b); }
function cdiv(a, b) { return Complex.from(a).div(b); }
function cexp(a) { return Complex.from(a).exp(); }
function csqrt(a) { return Complex.from(a).sqrt(); }

function mat2(A, B, Cc, D) {
  return [[Complex.from(A), Complex.from(B)], [Complex.from(Cc), Complex.from(D)]];
}

function matMul(A, B) {
  return [
    [A[0][0].mul(B[0][0]).add(A[0][1].mul(B[1][0])), A[0][0].mul(B[0][1]).add(A[0][1].mul(B[1][1]))],
    [A[1][0].mul(B[0][0]).add(A[1][1].mul(B[1][0])), A[1][0].mul(B[0][1]).add(A[1][1].mul(B[1][1]))]
  ];
}

function matVecMul(A, v) {
  return [
    A[0][0].mul(v[0]).add(A[0][1].mul(v[1])),
    A[1][0].mul(v[0]).add(A[1][1].mul(v[1]))
  ];
}

function matInv2(A) {
  const det = A[0][0].mul(A[1][1]).sub(A[0][1].mul(A[1][0]));
  return [
    [A[1][1].div(det), A[0][1].neg().div(det)],
    [A[1][0].neg().div(det), A[0][0].div(det)]
  ];
}

function linspace(start, stop, num) {
  if (num <= 1) return [start];
  const arr = [];
  const step = (stop - start) / (num - 1);
  for (let i = 0; i < num; i++) arr.push(start + step * i);
  return arr;
}

function arange(start, stop, step) {
  const arr = [];
  for (let x = start; x <= stop + step * 1e-9; x += step) arr.push(x);
  return arr;
}

function parseComplexToken(token) {
  // Accept examples: 3.42, 1, 400+500i, 400-500i, 2.1+0.03j
  const s = String(token).trim().replace(/j/g, "i").replace(/\s+/g, "");
  if (!s.includes("i")) return C(Number(s), 0);
  const noI = s.replace("i", "");
  const match = noI.match(/^([+-]?\d*\.?\d+(?:e[+-]?\d+)?)([+-]\d*\.?\d+(?:e[+-]?\d+)?)$/i);
  if (match) return C(Number(match[1]), Number(match[2]));
  if (noI === "" || noI === "+") return C(0, 1);
  if (noI === "-") return C(0, -1);
  return C(0, Number(noI));
}

function parseList(str) {
  return String(str)
    .split(/[;,\n]+/)
    .map(s => s.trim())
    .filter(Boolean);
}

function tmm1D(nList, dList, freqHz) {
  const k0 = 2 * Math.PI * freqHz / C0;
  const k = nList.map(n => Complex.from(n).mul(k0));

  const D0 = mat2(1, 1, k0, -k0);
  const D0inv = matInv2(D0);

  const D = [];
  const Dinv = [];
  const P = [];
  const Pinv = [];

  let Mat = null;

  for (let l = 0; l < k.length; l++) {
    D[l] = mat2(1, 1, k[l], k[l].neg());
    Dinv[l] = matInv2(D[l]);

    const phaseForward = I.mul(k[l]).mul(dList[l]).exp();
    const phaseBackward = I.neg().mul(k[l]).mul(dList[l]).exp();
    P[l] = mat2(phaseForward, 0, 0, phaseBackward);
    Pinv[l] = matInv2(P[l]);

    const segment = matMul(matMul(D[l], Pinv[l]), Dinv[l]);

    if (l === 0) {
      Mat = matMul(D0inv, segment);
    } else if (l === k.length - 1) {
      Mat = matMul(matMul(Mat, segment), D0);
    } else {
      Mat = matMul(Mat, segment);
    }
  }

  return { k, D, Dinv, P, Pinv, D0, Mat };
}

function calcTR(nList, dList, freqHz) {
  const { Mat } = tmm1D(nList, dList, freqHz);
  const t = C(1, 0).div(Mat[0][0]);
  const r = Mat[1][0].div(Mat[0][0]);
  return {
    T: t.abs() ** 2,
    R: r.abs() ** 2,
    t,
    r,
    Mat
  };
}

function tmm1DEfield(nList, dList, freqHz, pointsPerLayer = 500) {
  const { k, D, Dinv, P, D0, Mat } = tmm1D(nList, dList, freqHz);

  const vec0 = matVecMul(Mat, [C(1, 0), C(0, 0)]);
  const vec = [];
  vec[0] = matVecMul(matMul(Dinv[0], D0), vec0);

  let totalD = linspace(0.001e-6, dList[0], pointsPerLayer);
  let totalE = totalD.map(z => {
    const Ef = vec[0][0].mul(I.mul(k[0]).mul(z).exp());
    const Eb = vec[0][1].mul(I.neg().mul(k[0]).mul(z).exp());
    return Ef.add(Eb).abs();
  });

  for (let l = 1; l < nList.length; l++) {
    vec[l] = matVecMul(matMul(matMul(Dinv[l], D[l - 1]), P[l - 1]), vec[l - 1]);
    const zLocal = linspace(0, dList[l], pointsPerLayer);
    const offset = totalD[totalD.length - 1];
    const eLocal = zLocal.map(z => {
      const Ef = vec[l][0].mul(I.mul(k[l]).mul(z).exp());
      const Eb = vec[l][1].mul(I.neg().mul(k[l]).mul(z).exp());
      return Ef.add(Eb).abs();
    });
    totalE = totalE.concat(eLocal);
    totalD = totalD.concat(zLocal.map(z => z + offset));
  }

  return { totalD, totalE, height: Math.max(...totalE) };
}

function findPeaks(y, x) {
  const peaks = [];
  for (let i = 1; i < y.length - 1; i++) {
    if (y[i] > y[i - 1] && y[i] > y[i + 1]) {
      peaks.push({ x: x[i], y: y[i], index: i });
    }
  }
  return peaks;
}

function scanTMM({
  nList,
  dList,
  fStartTHz = 0,
  fStopTHz = 2.5,
  fStepTHz = 0.001,
  bandgapTHz = [0.3, 0.55],
  peakThresholdFraction = 0.001,
  peakNumber = 1
}) {
  const fTHz = arange(fStartTHz, fStopTHz, fStepTHz);
  const fHz = fTHz.map(v => v * 1e12);
  const T = [];
  const R = [];

  for (const f of fHz) {
    const out = calcTR(nList, dList, f);
    T.push(out.T);
    R.push(out.R);
  }

  const maxT = Math.max(...T);
  let peaks = findPeaks(T, fHz).filter(p => p.y > peakThresholdFraction * maxT);
  peaks = peaks.filter(p => p.x > bandgapTHz[0] * 1e12 && p.x < bandgapTHz[1] * 1e12);

  let field = null;
  if (peaks.length >= peakNumber) {
    field = tmm1DEfield(nList, dList, peaks[peakNumber - 1].x, 500);
  }

  return { fTHz, fHz, T, R, peaks, field };
}

function qwRefractiveIndex({ freqHz, B, mEff = 0.07 * ME, gamma = 5.9e9, nE = 10 * 3.2e11 * 1e4, dQW = 2.3e-6, epsBg = 3.6 ** 2 }) {
  const omega = 2 * Math.PI * freqHz;
  const omegaC = QE * B / mEff;
  const numerator = nE * QE ** 2 / mEff;
  const denom = C(gamma, -(omega - omegaC)); // gamma - i(omega - omega_c)
  const sigma = C(numerator, 0).div(denom);
  const epsCR = C(epsBg, 0).add(I.mul(sigma).div(EPS0 * omega * dQW));
  return epsCR.sqrt();
}

function defaultCavity() {
  const nAir = C(1, 0);
  const nSi = C(3.42, 0);
  const nDefect = C(3.42, 0);
  const dAir = 197e-6;
  const dSi = 50e-6;
  const dDefect = 100e-6;

  return {
    nList: [nSi, nAir, nSi, nAir, nDefect, nAir, nSi, nAir, nSi],
    dList: [dSi, dAir, dSi, dAir, dDefect, dAir, dSi, dAir, dSi]
  };
}

// Expose functions globally for simple HTML usage.
window.Complex = Complex;
window.C = C;
window.parseComplexToken = parseComplexToken;
window.parseList = parseList;
window.tmm1D = tmm1D;
window.calcTR = calcTR;
window.tmm1DEfield = tmm1DEfield;
window.scanTMM = scanTMM;
window.qwRefractiveIndex = qwRefractiveIndex;
window.defaultCavity = defaultCavity;
