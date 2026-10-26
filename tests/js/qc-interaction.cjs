/* Run the JavaScript embedded in a generated QC HTML file using only Node
 * built-ins. The DOM/canvas mock records interactions; it does not reproduce
 * the viewer's geometry or interaction implementation. Browser visual QA is
 * still required for actual painting and device-specific event delivery. */
"use strict";
const assert = require("assert").strict;
const fs = require("fs");
const vm = require("vm");

assert(process.argv[2], "Supply a generated QC HTML path");
const html = fs.readFileSync(process.argv[2], "utf8");
const sceneMatch = html.match(/<script\s+id="scene-data"\s+type="application\/json">([\s\S]*?)<\/script>/);
assert(sceneMatch, "Generated HTML contains its scene data");
const source = Array.from(html.matchAll(/<script([^>]*)>([\s\S]*?)<\/script>/g))
  .filter(match => !/application\/json/.test(match[1]))
  .map(match => match[2]).join("\n");
assert(source.trim(), "Generated HTML contains executable viewer JavaScript");

function close(actual, expected, message, tolerance = 1e-9) {
  assert(Number.isFinite(actual), `${message}: non-finite result ${actual}`);
  assert(Math.abs(actual - expected) <= tolerance,
    `${message}: expected ${expected}, received ${actual}`);
}
function vectorClose(actual, expected, message, tolerance = 1e-9) {
  assert.equal(actual.length, expected.length, `${message}: vector length`);
  actual.forEach((value, i) => close(value, expected[i], `${message}[${i}]`, tolerance));
}
const dot = (a, b) => a.reduce((sum, value, i) => sum + value * b[i], 0);
function sameOrientation(a, b, message) {
  // q and -q describe the same rotation.
  close(Math.abs(dot(a, b)), 1, message, 1e-10);
}
function changedOrientation(a, b, message) {
  assert(Math.abs(dot(a, b)) < 1 - 1e-8, message);
}

function makeViewer() {
  const scene = JSON.parse(sceneMatch[1]);
  const elements = new Map();
  const listeners = new Map();
  const captures = new Set();
  const frames = new Map();
  const counts = {draws: 0, rectangles: 0, captures: 0, releases: 0};
  const rectangle = {left: 37, top: 21, width: 800, height: 600,
    right: 837, bottom: 621};
  let nextFrame = 1;

  function element() {
    const attributes = new Map();
    const classes = new Set();
    return {
      textContent: "", innerHTML: "", style: {}, dataset: {}, children: [],
      hidden: false, width: 800, height: 600,
      classList: {add: x => classes.add(x), remove: x => classes.delete(x),
        contains: x => classes.has(x)},
      setAttribute: (name, value) => attributes.set(name, String(value)),
      getAttribute: name => attributes.get(name),
      appendChild(child) { this.children.push(child); return child; },
      querySelectorAll: () => []
    };
  }
  for (const match of html.matchAll(/\bid="([^"]+)"/g)) {
    elements.set(match[1], element());
  }
  elements.get("scene-data").textContent = JSON.stringify(scene);
  const canvas = elements.get("scene");
  const drawing = new Proxy({
    clearRect() { counts.draws++; },
    measureText(text) { return {width: String(text).length * 7}; }
  }, {get(target, property) {
    return property in target ? target[property] : () => {};
  }});
  canvas.getContext = () => drawing;
  canvas.getBoundingClientRect = () => { counts.rectangles++; return {...rectangle}; };
  canvas.setPointerCapture = id => { captures.add(id); counts.captures++; };
  canvas.hasPointerCapture = id => captures.has(id);
  canvas.releasePointerCapture = id => { captures.delete(id); counts.releases++; };
  canvas.addEventListener = (name, callback, options) => listeners.set(name, {callback, options});
  const requestAnimationFrame = callback => {
    const id = nextFrame++;
    frames.set(id, callback);
    return id;
  };
  const cancelAnimationFrame = id => frames.delete(id);
  const document = {
    getElementById(id) { assert(elements.has(id), `Unknown DOM element ${id}`); return elements.get(id); },
    querySelectorAll: () => [],
    createElement: () => element(),
    createTextNode: text => ({textContent: text})
  };
  const window = {devicePixelRatio: 1, requestAnimationFrame, cancelAnimationFrame};
  const context = vm.createContext({document, window, requestAnimationFrame,
    cancelAnimationFrame, console, performance: {now: () => 0}});
  vm.runInContext(source, context, {filename: "generated-qc-viewer.js"});
  const evaluate = code => vm.runInContext(code, context);
  const value = code => {
    const serialized = evaluate(`JSON.stringify(${code})`);
    assert.notEqual(serialized, undefined, `Generated viewer exposes ${code}`);
    return JSON.parse(serialized);
  };
  function frame() {
    const batch = Array.from(frames.values());
    frames.clear();
    batch.forEach(callback => callback(16));
  }
  function flush() {
    for (let i = 0; frames.size; i++) {
      assert(i < 10, "Viewer animation frames settle");
      frame();
    }
  }
  function pointer(type, id, x, y, extra = {}) {
    let prevented = false;
    const event = {pointerId: id, pointerType: "mouse", isPrimary: true,
      button: 0, buttons: type === "up" || type === "cancel" ? 0 : 1,
      clientX: rectangle.left + x, clientY: rectangle.top + y,
      preventDefault() { prevented = true; }, ...extra};
    const handler = canvas[`onpointer${type}`];
    assert.equal(typeof handler, "function", `pointer${type} handler exists`);
    handler(event);
    return prevented;
  }
  function drag(start, end, steps = 1, extra = {}) {
    pointer("down", 1, ...start, extra);
    for (let i = 1; i <= steps; i++) {
      pointer("move", 1, start[0] + (end[0] - start[0]) * i / steps,
        start[1] + (end[1] - start[1]) * i / steps, extra);
    }
    pointer("up", 1, ...end, extra);
  }
  function wheel(deltaY, deltaMode = 0) {
    const handler = listeners.get("wheel");
    assert(handler, "Wheel listener exists");
    assert.equal(handler.options.passive, false, "Wheel handler can prevent page scrolling");
    let prevented = false;
    handler.callback({deltaY, deltaMode, preventDefault() { prevented = true; }});
    assert(prevented, "Wheel gesture prevents page scrolling");
  }
  flush();
  return {evaluate, value, pointer, drag, wheel, frame, flush, frames, counts,
    captures, canvas, elements, rectangle, scene,
    orientation: () => value("state.orientation"),
    zoom: () => value("state.zoom"),
    project: point => value(`project(${JSON.stringify(point)})`)};
}

function checkRotation(viewer, message) {
  viewer.flush();
  const q = viewer.orientation();
  assert.equal(q.length, 4, `${message}: quaternion has four components`);
  close(Math.hypot(...q), 1, `${message}: normalized quaternion`, 1e-10);
  const center = viewer.scene.bounds.center;
  const origin = viewer.project(center);
  const scale = Math.min(viewer.rectangle.width, viewer.rectangle.height) /
    (viewer.scene.bounds.radius * 2.25) * viewer.zoom();
  // These RAS displacements correspond to the viewer's [x, z, -y] basis.
  const vectors = [[1, 0, 0], [0, 0, 1], [0, -1, 0]].map(delta => {
    const p = viewer.project(center.map((value, i) => value + delta[i]));
    return [(p.x - origin.x) / scale, -(p.y - origin.y) / scale, p.z - origin.z];
  });
  for (let i = 0; i < 3; i++) {
    for (let j = 0; j < 3; j++) close(dot(vectors[i], vectors[j]), i === j ? 1 : 0,
      `${message}: rotation basis ${i},${j}`, 1e-9);
  }
  const [a, b, c] = vectors;
  const determinant = a[0] * (b[1] * c[2] - b[2] * c[1]) -
    b[0] * (a[1] * c[2] - a[2] * c[1]) + c[0] * (a[1] * b[2] - a[2] * b[1]);
  close(determinant, 1, `${message}: no reflection or scaling`, 1e-9);
}

let passed = 0;
function test(name, check) {
  check();
  passed++;
  process.stdout.write(`PASS: ${name}\n`);
}

test("initial orientation preserves the previous RAS projection", () => {
  const viewer = makeViewer();
  const c = viewer.scene.bounds.center;
  const yaw = -.62, pitch = .24;
  const cy = Math.cos(yaw), sy = Math.sin(yaw), cp = Math.cos(pitch), sp = Math.sin(pitch);
  const scale = 600 / (viewer.scene.bounds.radius * 2.25);
  for (const delta of [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], [3, -4, 2]]) {
    const point = c.map((value, i) => value + delta[i]);
    const x = delta[0], y = delta[2], z = -delta[1];
    const x1 = cy * x + sy * z, z1 = -sy * x + cy * z;
    const y1 = cp * y - sp * z1, z2 = sp * y + cp * z1;
    const actual = viewer.project(point);
    vectorClose([actual.x, actual.y, actual.z], [400 + x1 * scale, 300 - y1 * scale, z2],
      `Previous projection at ${delta}`);
  }
  const expectedQuaternion = [Math.sin(pitch / 2) * Math.cos(yaw / 2),
    Math.cos(pitch / 2) * Math.sin(yaw / 2),
    Math.sin(pitch / 2) * Math.sin(yaw / 2),
    Math.cos(pitch / 2) * Math.cos(yaw / 2)];
  sameOrientation(viewer.orientation(), expectedQuaternion, "Initial Rx(pitch) Ry(yaw)");
  checkRotation(viewer, "Initial orientation");
});

test("vertical and edge drags keep rotating without pole clamps", () => {
  const viewer = makeViewer();
  let previous = viewer.orientation();
  for (let i = 0; i < 8; i++) {
    viewer.drag([400, 300], [400, 40], 5);
    changedOrientation(previous, viewer.orientation(), `Vertical drag ${i + 1} rotates`);
    checkRotation(viewer, `Vertical drag ${i + 1}`);
    previous = viewer.orientation();
  }
  viewer.drag([790, 10], [880, -80], 5);
  changedOrientation(previous, viewer.orientation(), "Captured drag beyond the canvas rotates");
  checkRotation(viewer, "Edge drag");
});

test("the same gesture is independent of pointer event frequency", () => {
  const one = makeViewer(), many = makeViewer();
  one.drag([230, 170], [670, 430], 1);
  many.drag([230, 170], [670, 430], 37);
  sameOrientation(one.orientation(), many.orientation(), "One versus 37 pointer moves");
  checkRotation(many, "High-frequency gesture");
});

test("untracked pointer IDs do not rotate or terminate the active gesture", () => {
  const viewer = makeViewer();
  const initial = viewer.orientation();
  viewer.pointer("down", 1, 300, 300);
  viewer.pointer("move", 99, 600, 100);
  viewer.pointer("up", 99, 600, 100);
  sameOrientation(viewer.orientation(), initial, "Unrelated move is ignored");
  assert(viewer.captures.has(1), "Unrelated pointer-up preserves active capture");
  viewer.pointer("move", 1, 440, 320);
  changedOrientation(initial, viewer.orientation(), "Original pointer still rotates");
  viewer.pointer("up", 1, 440, 320);
  assert(!viewer.captures.has(1), "Pointer-up releases capture");
});

test("cancel and lost capture end dragging and permit a fresh gesture", () => {
  for (const termination of ["cancel", "lostpointercapture"]) {
    const viewer = makeViewer();
    viewer.pointer("down", 1, 300, 300);
    viewer.pointer("move", 1, 380, 320);
    if (termination === "cancel") {
      viewer.pointer("cancel", 1, 380, 320);
      assert(!viewer.captures.has(1), "Cancel releases pointer capture");
    } else {
      viewer.captures.delete(1); // The browser has already revoked capture.
      assert.equal(typeof viewer.canvas.onlostpointercapture, "function");
      viewer.canvas.onlostpointercapture({pointerId: 1, pointerType: "mouse"});
    }
    assert(!viewer.canvas.classList.contains("dragging"), `${termination} clears dragging class`);
    const stopped = viewer.orientation();
    viewer.pointer("move", 1, 700, 500);
    sameOrientation(viewer.orientation(), stopped, `${termination} ignores later moves`);
    viewer.drag([300, 300], [450, 300]);
    changedOrientation(stopped, viewer.orientation(), `Fresh drag after ${termination}`);
  }
});

test("secondary mouse buttons cannot start rotation", () => {
  for (const button of [1, 2]) {
    const viewer = makeViewer();
    const initial = viewer.orientation();
    viewer.pointer("down", 1, 300, 300, {button, buttons: button === 1 ? 4 : 2});
    viewer.pointer("move", 1, 550, 450, {button, buttons: button === 1 ? 4 : 2});
    sameOrientation(viewer.orientation(), initial, `Mouse button ${button} ignored`);
    assert.equal(viewer.captures.size, 0, "Secondary mouse button does not capture");
  }
});

test("reset and resize cancel active gestures without a later pointer jump", () => {
  for (const action of ["reset", "resize"]) {
    const viewer = makeViewer();
    const initial = viewer.orientation();
    viewer.pointer("down", 1, 300, 300);
    viewer.pointer("move", 1, 440, 350);
    if (action === "reset") viewer.elements.get("reset").onclick();
    else viewer.evaluate("window.onresize()");
    assert.equal(viewer.captures.size, 0, `${action} releases active capture`);
    assert(!viewer.canvas.classList.contains("dragging"), `${action} clears dragging state`);
    const stable = viewer.orientation();
    if (action === "reset") sameOrientation(stable, initial, "Reset restores initial view immediately");
    viewer.pointer("move", 1, 700, 550);
    sameOrientation(viewer.orientation(), stable, `${action} ignores the interrupted gesture`);
    checkRotation(viewer, `After ${action} during a gesture`);
  }
});

test("two-touch pinch zoom transitions to one-touch rotation without a jump", () => {
  const viewer = makeViewer();
  const touch = {pointerType: "touch"};
  const orientation = viewer.orientation();
  viewer.pointer("down", 1, 300, 300, touch);
  viewer.pointer("down", 2, 500, 300, {...touch, isPrimary: false});
  viewer.pointer("move", 2, 600, 300, {...touch, isPrimary: false});
  close(viewer.zoom(), 1.5, "Pinch zoom follows distance ratio");
  sameOrientation(viewer.orientation(), orientation, "Pinch does not rotate");
  viewer.pointer("up", 2, 600, 300, {...touch, isPrimary: false});
  const afterPinch = viewer.orientation();
  viewer.pointer("move", 1, 300, 300, touch);
  sameOrientation(viewer.orientation(), afterPinch, "Remaining finger is rebased at its current position");
  close(viewer.zoom(), 1.5, "Lifting a finger preserves zoom");
  viewer.pointer("move", 1, 350, 320, touch);
  changedOrientation(afterPinch, viewer.orientation(), "Remaining finger resumes rotation");
  viewer.pointer("up", 1, 350, 320, touch);
  assert.equal(viewer.captures.size, 0, "Both touch captures released");
  checkRotation(viewer, "After pinch-to-rotation transition");
});

test("wheel delta modes are normalized and zoom remains bounded", () => {
  const pixels = makeViewer(), lines = makeViewer(), pages = makeViewer(), pagePixels = makeViewer();
  pixels.wheel(16, 0);
  lines.wheel(1, 1);
  close(lines.zoom(), pixels.zoom(), "One line equals 16 CSS pixels");
  pages.wheel(1, 2);
  pagePixels.wheel(600, 0);
  close(pages.zoom(), pagePixels.zoom(), "One page equals the canvas height");
  const viewer = makeViewer();
  viewer.wheel(1e6);
  close(viewer.zoom(), .35, "Minimum zoom");
  viewer.wheel(-1e6);
  close(viewer.zoom(), 5, "Maximum zoom");
});

test("bursts of renders and pointer events coalesce into one animation frame", () => {
  const viewer = makeViewer();
  viewer.counts.draws = 0;
  viewer.evaluate("for (let i = 0; i < 50; i++) render()");
  assert.equal(viewer.counts.draws, 0, "Render queues rather than drawing synchronously");
  assert.equal(viewer.frames.size, 1, "One pending frame for 50 render requests");
  viewer.frame();
  assert.equal(viewer.counts.draws, 1, "Exactly one draw for a burst");
  assert.equal(viewer.frames.size, 0, "No unnecessary continuous render loop");
  viewer.counts.draws = 0;
  viewer.drag([300, 300], [600, 380], 50);
  assert.equal(viewer.counts.draws, 0, "Pointer moves defer drawing");
  assert.equal(viewer.frames.size, 1, "One pending frame for pointer moves");
  viewer.frame();
  assert.equal(viewer.counts.draws, 1, "Pointer burst paints once");
});

test("canvas layout is measured per frame rather than per projected point", () => {
  const viewer = makeViewer();
  function drawWithPoints(n) {
    viewer.evaluate(`data.brain=Array.from({length:${n}},(_,i)=>[12+i%3,-7+i%5,5+i%7]);state.showLabels=true`);
    viewer.counts.rectangles = 0;
    viewer.evaluate("drawScene()");
    return viewer.counts.rectangles;
  }
  const small = drawWithPoints(12), large = drawWithPoints(12000);
  assert(small >= 1 && small <= 4, "Drawing reads a bounded number of canvas rectangles");
  assert.equal(large, small, "Adding points does not cause additional layout reads");
  viewer.counts.rectangles = 0;
  viewer.evaluate("for(let i=0;i<3000;i++) project([12,-7,5])");
  assert.equal(viewer.counts.rectangles, 0, "Projection uses cached frame geometry");
});

test("repeated rotations retain a proper rotation and reset restores the initial view", () => {
  const viewer = makeViewer();
  const initial = viewer.orientation();
  for (let i = 0; i < 240; i++) {
    const start = [180 + i % 5 * 30, 110 + i % 7 * 20];
    const end = [620 - i % 7 * 25, 490 - i % 5 * 25];
    viewer.drag(start, end, 3);
  }
  checkRotation(viewer, "After 240 gestures");
  viewer.wheel(-200);
  viewer.elements.get("reset").onclick();
  sameOrientation(viewer.orientation(), initial, "Reset restores initial orientation");
  close(viewer.zoom(), 1, "Reset restores zoom");
  checkRotation(viewer, "After reset");
});

process.stdout.write(`All generated QC interaction checks passed (${passed} groups).\n`);
