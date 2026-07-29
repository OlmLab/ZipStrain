// Minimal tab switcher for the parameter-selection page (self-contained; the readthedocs
// theme provides no tab component). Attaches to any .zs-tabs block.
(function () {
  function init() {
    document.querySelectorAll(".zs-tabs").forEach(function (box) {
      box.querySelectorAll(".zs-tab").forEach(function (btn) {
        btn.addEventListener("click", function () {
          box.querySelectorAll(".zs-tab").forEach(function (b) { b.classList.remove("active"); });
          box.querySelectorAll(".zs-panel").forEach(function (p) { p.classList.remove("active"); });
          btn.classList.add("active");
          var panel = box.querySelector("#" + btn.dataset.t);
          if (panel) panel.classList.add("active");
        });
      });
    });
  }
  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", init);
  } else {
    init();
  }
})();
