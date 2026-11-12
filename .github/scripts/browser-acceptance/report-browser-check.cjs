"use strict";
// Runs only on a native Windows runner. This is browser acceptance evidence,
// not anatomical validation or a substitute for inspecting the screenshots.
const assert = require("node:assert/strict");
const fs = require("node:fs");
const path = require("node:path");
const os = require("node:os");
const crypto = require("node:crypto");
const {pathToFileURL, fileURLToPath} = require("node:url");
const {chromium} = require("playwright");

assert(process.argv[2], "Supply browser-fixture.json from generate-fixture.R");
const manifestPath = path.resolve(process.argv[2]);
const manifest = JSON.parse(fs.readFileSync(manifestPath, "utf8"));
const root = fs.realpathSync(manifest.relocated_root);
const evidenceDir = path.join(path.dirname(manifestPath), "browser-evidence");
fs.mkdirSync(evidenceDir, {recursive: true});
const result = {
  passed: false,
  platform: process.platform,
  os_release: os.release(),
  package_version: manifest.package_version,
  source_commit: manifest.source_commit,
  fixture: manifest.fixture,
  checks: [], pages: [], downloads: [], screenshots: [], page_errors: [], image_loading: [],
  scope: "Synthetic reports in native Windows Chromium; manual and real-brain review remain separate."
};
let browser, context, page;
const hash = bytes => crypto.createHash("sha256").update(bytes).digest("hex");
function check(value, description, details) {
  result.checks.push({description, passed: Boolean(value), ...(details === undefined ? {} : {details})});
  assert(value, description);
}
function localFile(url) {
  const parsed = new URL(url);
  assert.equal(parsed.protocol, "file:", "Expected a local file URL");
  const file = fileURLToPath(parsed);
  const real = fs.realpathSync(file);
  const relative = path.relative(root, real);
  assert(relative && !relative.startsWith(`..${path.sep}`) && relative !== ".." && !path.isAbsolute(relative),
    `Local link must remain inside the relocated bundle: ${url}`);
  assert(fs.statSync(real).isFile(), `Local link must identify a file: ${url}`);
  return real;
}
const reportURL = kind => pathToFileURL(path.join(root, manifest.reports[kind])).href;
async function settled() {
  await page.evaluate(async () => {
    if (document.fonts) await document.fonts.ready;
    await new Promise(resolve => requestAnimationFrame(() => requestAnimationFrame(resolve)));
  });
}
async function open(url) {
  await page.goto(url, {waitUntil: "load"});
  await settled();
}
async function screenshot(name, fullPage = false) {
  const filename = `${name}.png`;
  await page.screenshot({path: path.join(evidenceDir, filename), fullPage});
  result.screenshots.push(filename);
}
async function exposeReportImages(relative) {
  const initiallyUnloaded = await page.locator("img").evaluateAll(images => images
    .filter(image => !image.complete || image.naturalWidth === 0)
    .map(image => ({url: image.src, complete: image.complete,
      inside_closed_details: Boolean(image.closest("details:not([open])"))})));
  // Chromium can defer images inside collapsed details. Exercise the same
  // disclosure and scrolling actions a reader uses before requiring decoding.
  // The assertions below still reject missing, corrupt or incorrectly linked images.
  const disclosures = page.locator("details");
  let expanded = 0;
  for (let i = 0; i < await disclosures.count(); i++) {
    const section = disclosures.nth(i);
    if (await section.locator("img").count() > 0 &&
        !await section.evaluate(element => element.open)) {
      await section.locator(":scope > summary").click();
      expanded++;
    }
  }
  const images = page.locator("img");
  for (let i = 0; i < await images.count(); i++) {
    await images.nth(i).scrollIntoViewIfNeeded();
  }
  await page.waitForFunction(() => [...document.images].every(image => image.complete));
  const decoded = await images.evaluateAll(elements => elements.map(image => ({
    url: image.src, complete: image.complete,
    natural_width: image.naturalWidth, natural_height: image.naturalHeight
  })));
  result.image_loading.push({relative_path: relative, initially_unloaded: initiallyUnloaded,
    expanded_sections: expanded, decoded});
  await page.evaluate(() => window.scrollTo(0, 0));
  await settled();
}
async function observeFileClick(link, targetURL, allowImageNavigation) {
  const startingURL = page.url();
  let finish;
  // Resolve every observation, including timeouts, so a failed click cannot
  // leave a rejecting event-wait promise behind. Remove both listeners once
  // either supported browser behavior has been observed.
  const observation = new Promise(resolve => {
    let completed = false, timer;
    const onDownload = download => finish({behavior: "download", download});
    const onNavigation = frame => {
      if (frame !== page.mainFrame() || frame.url() === startingURL) return;
      if (allowImageNavigation && frame.url() === targetURL) {
        finish({behavior: "image_navigation"});
      } else {
        finish({error: new Error(`Unexpected navigation after file click: ${frame.url()}`)});
      }
    };
    finish = outcome => {
      if (completed) return;
      completed = true;
      clearTimeout(timer);
      page.off("download", onDownload);
      page.off("framenavigated", onNavigation);
      resolve(outcome);
    };
    page.on("download", onDownload);
    page.on("framenavigated", onNavigation);
    timer = setTimeout(() => finish({error: new Error(
      `No supported file-open behavior observed within 15000ms: ${targetURL}`
    )}), 15000);
  });
  try {
    await link.click();
  } catch (error) {
    finish({error});
  }
  const outcome = await observation;
  if (outcome.error) throw outcome.error;
  return outcome;
}
async function canvasPixels() {
  return page.locator("#scene").evaluate(canvas => {
    const pixels = canvas.getContext("2d").getImageData(0, 0, canvas.width, canvas.height).data;
    let centralPainted = 0, centralFingerprint = 2166136261;
    const centralColors = new Set();
    for (let y = Math.floor(canvas.height * .2); y < canvas.height * .8; y++) {
      for (let x = Math.floor(canvas.width * .2); x < canvas.width * .8; x++) {
        const i = (y * canvas.width + x) * 4;
        for (let channel = 0; channel < 4; channel++) {
          centralFingerprint = Math.imul(centralFingerprint ^ pixels[i + channel], 16777619) >>> 0;
        }
        if (pixels[i + 3] > 0) {
          centralPainted++;
          centralColors.add(`${pixels[i]},${pixels[i + 1]},${pixels[i + 2]},${pixels[i + 3]}`);
        }
      }
    }
    return {width: canvas.width, height: canvas.height, centralPainted,
      centralColors: centralColors.size, centralFingerprint, image: canvas.toDataURL("image/png")};
  });
}
async function projectedSceneBounds() {
  // Inspect the coordinates actually projected by the generated viewer,
  // independently of whether its canvas happens to contain some painted pixels.
  return page.evaluate(() => {
    const points = [...data.brain, ...data.focal];
    for (const target of data.targets) {
      points.push(target.point, ...(target.boundary || []));
    }
    for (const streamline of data.streamlines) {
      points.push(...streamline.points);
      for (const endpoint of streamline.endpoints || []) points.push(endpoint.point);
    }
    const projected = points.map(point => project(point));
    const rectangle = canvas.getBoundingClientRect();
    return {point_count: projected.length, canvas_width: rectangle.width, canvas_height: rectangle.height,
      finite: projected.every(point => Number.isFinite(point.x) && Number.isFinite(point.y)),
      min_x: Math.min(...projected.map(point => point.x)), max_x: Math.max(...projected.map(point => point.x)),
      min_y: Math.min(...projected.map(point => point.y)), max_y: Math.max(...projected.map(point => point.y))};
  });
}

(async () => {
  try {
    check(process.platform === "win32", "Acceptance browser runs on native Windows");
    check(manifest.os === "Windows", "Synthetic reports were generated on Windows");
    check(manifest.copy_checksums_identical && manifest.original_outputs_preserved,
      "Relocated bundle preserves the generated original and identical file bytes");
    check(manifest.scientific_tables_passed === 7 && manifest.declared_values_passed === 18,
      "Installed-package fixture passed all scientific known answers");
    browser = await chromium.launch({headless: true});
    result.browser_version = browser.version();
    result.playwright_version = require("playwright/package.json").version;
    context = await browser.newContext({viewport: {width: 1440, height: 1000},
      deviceScaleFactor: 1, acceptDownloads: true});
    context.setDefaultTimeout(15000);
    page = await context.newPage();
    page.on("pageerror", error => result.page_errors.push(error.message));
    // Reports are self-contained. Do not contact external sites during review.
    await context.route(/^https?:/, route => route.abort());

    const localLinks = [];
    const ids = new Map();
    for (const [i, relative] of manifest.html_files.entries()) {
      const url = pathToFileURL(path.join(root, relative)).href;
      await open(url);
      await exposeReportImages(relative);
      const content = await page.evaluate(() => ({
        title: document.title,
        ids: [...document.querySelectorAll("[id]")].map(element => element.id),
        names: [...document.querySelectorAll("a[name]")].map(element => element.name),
        links: [...document.querySelectorAll("a[href]")].map(element => ({raw: element.getAttribute("href"), url: element.href})),
        images: [...document.images].map(element => ({url: element.src,
          loaded: element.complete && element.naturalWidth > 0 && element.naturalHeight > 0}))
      }));
      check(Boolean(content.title.trim()), `HTML title is present: ${relative}`);
      check(content.images.every(image => image.loaded), `Every report image loaded: ${relative}`,
        content.images.filter(image => !image.loaded));
      ids.set(url, new Set([...content.ids, ...content.names]));
      for (const link of content.links) {
        check(!link.raw.includes("\\"), `HTML link uses URL separators: ${relative}`);
        if (new URL(link.url).protocol === "file:") {
          localFile(link.url);
          localLinks.push({source: relative, url: link.url});
        }
      }
      for (const image of content.images) {
        if (new URL(image.url).protocol === "file:") localFile(image.url);
      }
      result.pages.push({relative_path: relative, image_count: content.images.length,
        local_link_count: content.links.filter(link => new URL(link.url).protocol === "file:").length});
      // Retain a viewport for long guides; the four primary views also get
      // complete screenshots below. Existence alone is not visual approval.
      await screenshot(`page-${String(i + 1).padStart(2, "0")}`);
      if (relative === manifest.reports.results) await screenshot("results-expanded", true);
    }
    for (const link of localLinks) {
      const target = new URL(link.url);
      const fragment = decodeURIComponent(target.hash.slice(1));
      target.hash = "";
      if (fragment && /[.]html$/i.test(target.pathname)) {
        check(ids.has(target.href) && ids.get(target.href).has(fragment),
          `Local HTML fragment resolves: ${link.source} -> ${fragment}`);
      }
    }
    result.local_references_checked = localLinks.length;
    check(localLinks.length > 0, "Report bundle contains working local navigation");

    await open(reportURL("index"));
    await screenshot("index", true);
    await page.getByRole("link", {name: "Results", exact: true}).click();
    await page.waitForURL(reportURL("results"));
    await settled();
    check(await page.locator("img").count() > 0, "Results report contains rendered figures");
    await screenshot("results", true);
    for (const kind of ["target", "viewer"]) {
      await open(reportURL("index"));
      const sections = page.locator("details");
      for (let i = 0; i < await sections.count(); i++) {
        if (!await sections.nth(i).evaluate(element => element.open)) {
          await sections.nth(i).locator("summary").first().click();
        }
      }
      const links = page.locator("a[href]");
      const position = await links.evaluateAll((elements, wanted) => elements.findIndex(element => element.href === wanted), reportURL(kind));
      check(position >= 0, `Index provides the ${kind} report link`);
      await links.nth(position).click();
      await page.waitForURL(reportURL(kind));
      await settled();
      await screenshot(kind, true);
    }

    await open(reportURL("target"));
    for (const extension of ["csv", "png"]) {
      const link = page.locator(`a[download][href$=".${extension}"]`).first();
      check(await link.count() === 1, `Target report offers a ${extension.toUpperCase()} download`);
      const url = await link.evaluate(element => element.href);
      const expected = fs.readFileSync(localFile(url));
      const observed = await observeFileClick(link, url, extension === "png");
      if (observed.behavior === "download") {
        const download = observed.download;
        check(await download.failure() === null, `${extension.toUpperCase()} download completes`);
        const received = fs.readFileSync(await download.path());
        check(hash(received) === hash(expected), `${extension.toUpperCase()} download matches the linked file`);
        result.downloads.push({extension, behavior: "download", target_url: url,
          page_url: page.url(), filename: download.suggestedFilename(),
          target_sha256: hash(expected), sha256: hash(received)});
      } else {
        // Chromium can display a local PNG despite an HTML download attribute.
        // Accept only the observed exact-target navigation and a decoded image
        // matching the dimensions declared by that linked PNG's IHDR header.
        check(extension === "png", "Only PNG links may open as browser images");
        check(expected.length >= 24 && expected.subarray(0, 8).equals(Buffer.from([137, 80, 78, 71, 13, 10, 26, 10])) &&
          expected.toString("ascii", 12, 16) === "IHDR", "Linked image has a valid PNG signature and dimension header");
        const expectedWidth = expected.readUInt32BE(16), expectedHeight = expected.readUInt32BE(20);
        await page.waitForLoadState("load");
        await page.waitForFunction(() => document.images.length === 1 && document.images[0].complete &&
          document.images[0].naturalWidth > 0 && document.images[0].naturalHeight > 0);
        const image = await page.evaluate(() => ({url: location.href, content_type: document.contentType,
          natural_width: document.images[0].naturalWidth, natural_height: document.images[0].naturalHeight}));
        check(image.url === url && image.content_type === "image/png", "PNG opens at the exact linked URL with image/png content type");
        check(image.natural_width === expectedWidth && image.natural_height === expectedHeight,
          "Browser-decoded PNG dimensions match the linked file header");
        await screenshot("png-opened-in-browser", true);
        result.downloads.push({extension, behavior: "image_navigation", target_url: url,
          target_sha256: hash(expected), page_url: image.url, content_type: image.content_type,
          natural_width: image.natural_width, natural_height: image.natural_height});
      }
    }

    await open(reportURL("viewer"));
    const initial = await canvasPixels();
    check(initial.width > 0 && initial.height > 0 && initial.centralPainted > 20 && initial.centralColors > 8,
      "3-D canvas paints nonuniform scene content away from the compass",
      {central_painted_pixels: initial.centralPainted, central_colors: initial.centralColors});
    await screenshot("viewer-initial");
    const bounds = await page.locator("#scene").boundingBox();
    const x = bounds.x + bounds.width * .55, y = bounds.y + bounds.height * .45;
    await page.mouse.move(x, y);
    await page.mouse.down();
    await page.mouse.move(x + 120, y + 70, {steps: 12});
    await page.mouse.up();
    await settled();
    const rotated = await canvasPixels();
    check(rotated.centralFingerprint !== initial.centralFingerprint,
      "Dragging changes central scene pixels, not only the compass");
    await screenshot("viewer-rotated");
    await page.mouse.wheel(0, -300);
    await settled();
    const zoomed = await canvasPixels();
    check(zoomed.centralFingerprint !== rotated.centralFingerprint,
      "Mouse wheel changes rendered zoom in the scene centre");
    await screenshot("viewer-zoomed");
    result.canvas = {width: initial.width, height: initial.height,
      central_painted_pixels: initial.centralPainted, central_colors: initial.centralColors,
      initial_fingerprint: initial.centralFingerprint, rotated_fingerprint: rotated.centralFingerprint,
      zoomed_fingerprint: zoomed.centralFingerprint};
    await page.locator("#reset").click();
    await settled();
    check((await canvasPixels()).image === initial.image, "Reset restores the initial rendered scene exactly");
    await screenshot("viewer-reset");

    // Drag near the right edge to exercise trackball roll with real pointer events.
    const edge = await page.locator("#scene").boundingBox();
    await page.mouse.move(edge.x + edge.width - 25, edge.y + edge.height * .45);
    await page.mouse.down();
    await page.mouse.move(edge.x + edge.width - 25, edge.y + edge.height * .65, {steps: 12});
    await page.mouse.up();
    await settled();
    check((await canvasPixels()).centralFingerprint !== initial.centralFingerprint,
      "Edge dragging rolls the rendered scene");
    await screenshot("viewer-rolled");
    await page.locator("#reset").click();
    await settled();
    check((await canvasPixels()).image === initial.image, "Reset restores the scene after edge roll");

    // Hide endpoint markers so they cannot completely occlude the tiny
    // synthetic focal boundary while testing that layer's visibility control.
    await page.locator("#toggleEndpoints").click();
    await settled();
    const layerBaseline = await canvasPixels();
    await page.locator("#toggleFocal").click();
    await settled();
    check(await page.locator("#toggleFocal").getAttribute("aria-pressed") === "false",
      "Focal-boundary toggle exposes its disabled state");
    check((await canvasPixels()).centralFingerprint !== layerBaseline.centralFingerprint,
      "Focal-boundary toggle changes central rendered pixels");
    await page.locator("#toggleFocal").click();
    await page.locator("#toggleEndpoints").click();
    await page.locator("#advancedControls > summary").click();
    await page.locator('.segment[data-mode="focal"]').click();
    await settled();
    check(await page.locator('.segment[data-mode="focal"]').getAttribute("aria-pressed") === "true",
      "Advanced path display selects the focal segment");
    check((await canvasPixels()).centralFingerprint !== initial.centralFingerprint,
      "Focal-segment selection changes the rendered paths");
    await page.locator('.segment[data-mode="points"]').click();
    await page.locator('.endpoint-mode[data-endpoint-mode="Dilated1"]').click();
    await settled();
    check(await page.locator('.endpoint-mode[data-endpoint-mode="Dilated1"]').getAttribute("aria-pressed") === "true",
      "Endpoint sensitivity control changes its selected state");
    check((await page.locator("#endpointStatus").innerText()).includes("1-voxel"),
      "Endpoint status reports the selected sensitivity mode");
    await screenshot("viewer-advanced");
    await page.locator('.endpoint-mode[data-endpoint-mode="Exact"]').click();
    await page.locator("#advancedControls > summary").click();
    await page.setViewportSize({width: 500, height: 900});
    await settled();
    await page.locator("#reset").click();
    await settled();
    check(await page.evaluate(() => document.documentElement.scrollWidth <= window.innerWidth + 2),
      "Narrow viewer fits the viewport without horizontal overflow");
    const narrow = await canvasPixels();
    check(narrow.centralPainted > 20 && narrow.centralColors > 8, "Narrow viewer still paints nonuniform scene content");
    const narrowBounds = await projectedSceneBounds();
    result.narrow_scene_bounds = narrowBounds;
    const markerPadding = 9; // Largest default marker: radius 7 plus half its 3-pixel stroke.
    check(narrowBounds.point_count > 0 && narrowBounds.finite &&
      narrowBounds.min_x >= markerPadding && narrowBounds.min_y >= markerPadding &&
      narrowBounds.max_x <= narrowBounds.canvas_width - markerPadding &&
      narrowBounds.max_y <= narrowBounds.canvas_height - markerPadding,
    "After reset, the complete projected scene fits inside the narrow canvas with marker padding", narrowBounds);
    await screenshot("viewer-narrow", true);
    check(result.page_errors.length === 0, "No report or viewer JavaScript errors", result.page_errors);
    result.passed = true;
  } catch (error) {
    result.error = error.stack || String(error);
    if (page) {
      result.failure_url = page.url();
      try { await screenshot("failure"); } catch (_) { /* Preserve the primary failure. */ }
    }
    process.exitCode = 1;
  } finally {
    fs.writeFileSync(path.join(evidenceDir, "browser-check.json"), JSON.stringify(result, null, 2) + "\n");
    if (browser) await browser.close();
    process.stdout.write(JSON.stringify({passed: result.passed, checks: result.checks.length,
      pages: result.pages.length, local_references: result.local_references_checked,
      error: result.error || null}, null, 2) + "\n");
  }
})();
