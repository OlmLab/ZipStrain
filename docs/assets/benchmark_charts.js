document.addEventListener("DOMContentLoaded", () => {
  const charts = [...document.querySelectorAll("img[data-release-chart]")];
  if (!charts.length) return;

  // Keep the bundled baseline when browsing offline or before the first release.
  Promise.all(charts.map((chart) => new Promise((resolve, reject) => {
    const releaseChart = new Image();
    releaseChart.onload = () => resolve(releaseChart.src);
    releaseChart.onerror = reject;
    releaseChart.src = chart.dataset.releaseChart;
  }))).then((sources) => {
    charts.forEach((chart, index) => { chart.src = sources[index]; });
    document.querySelector("#benchmark-source").textContent =
      "Release history from the GitHub Actions benchmark runner.";
  }).catch(() => {});
});
