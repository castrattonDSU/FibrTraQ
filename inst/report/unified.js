(function () {
  "use strict";
  const select = document.getElementById("ftq-run");
  const active = document.getElementById("ftq-active");
  if (!select || !active) return;

  // Run zero is already rendered in the HTML before the browser restores any
  // saved form selection on reload or return navigation.
  let current = "0";
  let request = 0;

  function resetViewer() {
    request++;
    const host = document.getElementById("ftq-viewer-host");
    if (host) host.replaceChildren();
    const button = document.getElementById("ftq-open-viewer");
    if (button) {
      button.disabled = false;
      button.hidden = false;
    }
  }

  function showSelectedRun() {
    if (select.value === current) return;
    resetViewer();
    const previous = document.getElementById("ftq-run-" + current);
    const next = document.getElementById("ftq-run-" + select.value);
    // Inactive analyses are inert. Only the selected analysis loads images.
    previous.content.append(...active.childNodes);
    active.replaceChildren(next.content);
    current = select.value;
  }
  select.addEventListener("change", showSelectedRun);
  window.addEventListener("pageshow", showSelectedRun);
  showSelectedRun();

  active.addEventListener("click", function (event) {
    const button = event.target.closest("#ftq-open-viewer");
    if (!button || button.disabled) return;

    const host = document.getElementById("ftq-viewer-host");
    const token = ++request;
    button.disabled = true;
    host.textContent = "Loading 3D scene…";
    const resource = document.createElement("script");
    // A classic local script works when a relocated report opens via file://.
    // fetch() of a neighbouring local JSON file is blocked by many browsers.
    resource.src = button.dataset.scene;
    function showError() {
      if (token !== request) return;
      host.textContent = "The 3D scene could not be loaded. Keep the figures folder with this report and try again.";
      button.disabled = false;
    }
    resource.onerror = function () {
      resource.remove();
      showError();
    };
    resource.onload = function () {
      const scene = resource.ftqScene;
      resource.remove();
      if (token !== request) return;
      if (!scene || !scene.metadata || !Array.isArray(scene.streamlines)) {
        showError();
        return;
      }
      // Re-encode '<' after serializing the JSON for the iframe's script tag.
      const json = JSON.stringify(scene).replace(/</g, "\\u003c");
      const frame = document.createElement("iframe");
      frame.id = "ftq-viewer";
      frame.title = scene.metadata.title;
      frame.srcdoc = document.getElementById("ftq-viewer-shell").content.textContent
        .replace("FTQ_EMBEDDED_SCENE_JSON", function () { return json; });
      host.replaceChildren(frame);
      button.hidden = true;
    };
    document.head.appendChild(resource);
  });
})();
