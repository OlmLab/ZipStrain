// Recombination-mosaics guide: size the embedded explorer and let links set its sliders.
(function () {
  function init() {
    var frame = document.getElementById("mosaic-explorer");
    if (!frame) return;
    var base = frame.getAttribute("src").split("?")[0];

    // The explorer reports its height, so the frame grows with it instead of scrolling.
    window.addEventListener("message", function (event) {
      var data = event.data;
      if (event.source === frame.contentWindow && data && data.type === "mosaic-explorer-height" && data.height > 200) {
        frame.style.height = (data.height + 8) + "px";
      }
    });

    // Links such as interactive/mosaic_explorer.html#factor=7&switch=0.13 drive the embedded copy.
    var clicks = 0;
    var links = document.querySelectorAll('a[href^="' + base + '#"]');
    Array.prototype.forEach.call(links, function (link) {
      link.addEventListener("click", function (event) {
        event.preventDefault();
        // A changing counter makes repeated clicks on one link re-apply its settings.
        clicks += 1;
        frame.src = base + "?embed=1" + link.hash + "&n=" + clicks;
        frame.parentElement.scrollIntoView({ behavior: "smooth", block: "start" });
      });
    });
  }
  if (document.readyState === "loading") document.addEventListener("DOMContentLoaded", init);
  else init();
})();
