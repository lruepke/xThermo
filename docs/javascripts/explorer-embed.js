/* Homepage: the xThermo Explorer (app/index.html?embed) in an iframe that takes the height of its content, so the
   page scrolls as one, and that follows the documentation's light/dark scheme. Same origin; no messages needed. */
(function () {
  function init() {
    var frame = document.getElementById("xt-explorer-frame");
    if (!frame) return;
    var scheme = function () {
      return document.body.getAttribute("data-md-color-scheme") === "slate" ? "dark" : "light";
    };
    var observer = null;
    function attach() {
      var doc;
      try { doc = frame.contentDocument; } catch (e) { return; }
      if (!doc || !doc.body) return;
      doc.documentElement.setAttribute("data-theme", scheme());
      var fit = function () { frame.style.height = Math.ceil(doc.body.getBoundingClientRect().height) + "px"; };
      if (observer) observer.disconnect();
      observer = new ResizeObserver(fit);
      observer.observe(doc.body);
      fit();
    }
    frame.addEventListener("load", attach);
    if (frame.contentDocument && frame.contentDocument.readyState === "complete") attach();
    new MutationObserver(function () {
      try { frame.contentDocument.documentElement.setAttribute("data-theme", scheme()); } catch (e) { /* ignore */ }
    }).observe(document.body, { attributes: true, attributeFilter: ["data-md-color-scheme"] });
  }
  if (document.readyState === "loading") document.addEventListener("DOMContentLoaded", init);
  else init();
})();
