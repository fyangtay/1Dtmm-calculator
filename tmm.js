function runTMM() {
  const lambda = Number(document.getElementById("lambda").value);

  // Temporary placeholder calculation
  const fakeReflectance = Math.sin(lambda / 100) ** 2;

  document.getElementById("result").innerText =
    "Example reflectance R = " + fakeReflectance.toFixed(4);
}
