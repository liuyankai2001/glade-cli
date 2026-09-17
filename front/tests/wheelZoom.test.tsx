import { afterEach, expect, it, vi } from "vitest";
import { act, cleanup, fireEvent, render, screen } from "@testing-library/react";
import "@testing-library/jest-dom/vitest";
import "../src/style.css";
import "../src/controls.css";
import { PlasmidRing } from "../src/PlasmidRing";
import type { Preview } from "../src/types";

const preview: Preview = {
  valid: true, issues: [], warnings: [], sequence: "", length_bp: 1000,
  gc_percent: 50, sequence_sha256: "sequence", source_fingerprint: "source",
  manifest_revision: 1, enzymes: [],
  segments: [
    { id: "amp", label: "AmpR", kind: "resistance", component_type: "resistance", instance_id: "amp", start_bp: 1, end_bp: 400, length_bp: 400 },
    { id: "ori", label: "ori", kind: "replication", component_type: "replication", instance_id: "ori", start_bp: 401, end_bp: 600, length_bp: 200 },
    { id: "expression", label: "expression", kind: "expression", component_type: "expression", instance_id: "expression", start_bp: 601, end_bp: 1000, length_bp: 400 },
  ],
  features: [{ label: "reverse gene", type: "CDS", kind: "gene", start_bp: 650, end_bp: 850, strand: -1 }],
};

function setup() {
  const reorder = vi.fn();
  const mounted = render(<section className="canvas">
    <PlasmidRing preview={preview} construct={null} enzymes={[]}
      componentOrder={["amp", "ori", "expression"]} onReorder={reorder} />
  </section>);
  const svg = screen.getByRole("img", { name: "质粒环图" });
  vi.spyOn(svg, "getBoundingClientRect").mockReturnValue({
    left: 0, top: 0, width: 440, height: 440,
    right: 440, bottom: 440, x: 0, y: 0, toJSON() {},
  } as DOMRect);
  const part = mounted.container.querySelector('[data-ring-instance="amp"]')!;
  const layer = () => part.parentElement!;
  const zoom = () => Number(layer().getAttribute("transform")!.match(/scale\(([^)]+)\)/)![1]);
  return { ...mounted, svg, part, layer, zoom, reorder };
}

function wheel(target: Element, deltaY: number, extra: WheelEventInit = {}) {
  const event = new WheelEvent("wheel", { bubbles: true, cancelable: true, deltaY, ...extra });
  fireEvent(target, event);
  return event;
}

afterEach(() => { cleanup(); vi.restoreAllMocks(); });

it("zooms in and out by wheel distance without changing DNA coordinates or order", () => {
  const { svg, part, zoom, reorder } = setup();
  expect(wheel(svg, -100).defaultPrevented).toBe(true);
  expect(zoom()).toBeCloseTo(1.1);
  wheel(svg, 100);
  expect(zoom()).toBeCloseTo(1);
  wheel(svg, -50);
  expect(zoom()).toBeCloseTo(1.048808848);
  expect(part.querySelector("path")).toHaveAttribute("data-start-bp", "1");
  expect(part.querySelector("path")).toHaveAttribute("data-end-bp", "400");
  expect(reorder).not.toHaveBeenCalled();
});

it("clamps at 50% and 200% and scales the central mask and bp number with the ring", () => {
  const { svg, zoom, layer, container } = setup();
  wheel(svg, -10000);
  expect(zoom()).toBe(2);
  wheel(svg, -100);
  expect(zoom()).toBe(2);
  wheel(svg, 10000);
  expect(zoom()).toBe(0.5);
  wheel(svg, 100);
  expect(zoom()).toBe(0.5);
  for (const node of [container.querySelector('circle[r="78"]'), container.querySelector(".ring-number"), container.querySelector(".ring-unit")]) {
    expect(node!.closest("g[transform]")).toBe(layer());
  }
});

it("normalizes line and page wheel units to the same pixel distance", () => {
  const { svg, zoom } = setup();
  Object.defineProperty(svg, "clientHeight", { configurable: true, value: 400 });
  wheel(svg, -80);
  const pixels = zoom();
  expect(pixels).toBeCloseTo(1.079230345);
  fireEvent.click(screen.getByRole("button", { name: "重置环图大小" }));
  wheel(svg, -5, { deltaMode: 1 });
  expect(zoom()).toBeCloseTo(pixels);
  fireEvent.click(screen.getByRole("button", { name: "重置环图大小" }));
  wheel(svg, -0.2, { deltaMode: 2 });
  expect(zoom()).toBeCloseTo(pixels);
});

it.each([{ ctrlKey: true }, { metaKey: true }, { deltaX: 100 }])("keeps browser wheel behavior for %j", (extra) => {
  const { svg, zoom } = setup();
  const event = wheel(svg, "deltaX" in extra ? 0 : -100, extra);
  expect(event.defaultPrevented).toBe(false);
  expect(zoom()).toBe(1);
});

it("leaves scrolling outside the SVG and over its buttons alone", () => {
  const { zoom, container } = setup();
  expect(wheel(container.querySelector(".canvas")!, -100).defaultPrevented).toBe(false);
  expect(wheel(screen.getByRole("button", { name: "显示基因注释" }), -100).defaultPrevented).toBe(false);
  expect(zoom()).toBe(1);
});

it("resets to 100% while retaining module direction and gene annotation state", () => {
  const { svg, part, zoom } = setup();
  fireEvent.click(part);
  fireEvent.click(screen.getByRole("button", { name: "显示基因注释" }));
  wheel(svg, -10000);
  expect(zoom()).toBe(2);
  fireEvent.click(screen.getByRole("button", { name: "重置环图大小" }));
  expect(zoom()).toBe(1);
  expect(part).toHaveAttribute("aria-pressed", "true");
  expect(screen.getByTestId("module-5-prime")).toHaveAttribute("data-bp", "1");
  expect(screen.getByTestId("module-3-prime")).toHaveAttribute("data-bp", "400");
  expect(screen.getByTestId("strand-arrow")).toHaveAttribute("data-strand", "-1");
  expect(screen.getByRole("button", { name: "隐藏基因注释" })).toHaveAttribute("aria-pressed", "true");
});

it("cancels a live drag on zoom and suppresses its release click", () => {
  const { svg, part, reorder } = setup();
  fireEvent.mouseDown(part, { button: 0, clientX: 220, clientY: 80 });
  fireEvent.mouseMove(window, { clientX: 250, clientY: 80 });
  expect(screen.getByTestId("ring-drag-ghost")).toBeInTheDocument();
  wheel(svg, -100);
  expect(screen.queryByTestId("ring-drag-ghost")).not.toBeInTheDocument();
  fireEvent.mouseUp(part, { clientX: 250, clientY: 80 });
  fireEvent.click(part, { detail: 1, clientX: 250, clientY: 80 });
  expect(reorder).not.toHaveBeenCalled();
  expect(screen.queryByTestId("module-direction")).not.toBeInTheDocument();
});

it.each(["wheel", "reset"])("cancels dragging before an immediate mouseup after %s", (action) => {
  const { svg, part, reorder, zoom } = setup();
  if (action === "reset") wheel(svg, -100);
  fireEvent.mouseDown(part, { button: 0, clientX: 220, clientY: 80 });
  fireEvent.mouseMove(window, { clientX: 250, clientY: 80 });
  expect(screen.getByTestId("ring-drag-ghost")).toBeInTheDocument();
  act(() => {
    if (action === "wheel") svg.dispatchEvent(new WheelEvent("wheel", { bubbles: true, cancelable: true, deltaY: -100 }));
    else screen.getByRole("button", { name: "重置环图大小" }).dispatchEvent(new MouseEvent("click", { bubbles: true }));
    window.dispatchEvent(new MouseEvent("mouseup", { clientX: 250, clientY: 80 }));
  });
  expect(reorder).not.toHaveBeenCalled();
  expect(zoom()).toBeCloseTo(action === "wheel" ? 1.1 : 1);
  expect(screen.queryByTestId("ring-drag-ghost")).not.toBeInTheDocument();
  fireEvent.click(part, { detail: 1, clientX: 250, clientY: 80 });
  expect(screen.queryByTestId("module-direction")).not.toBeInTheDocument();
});

it("composes wheel events received before a render", () => {
  const { svg, zoom } = setup();
  act(() => {
    for (let i = 0; i < 3; i++) svg.dispatchEvent(new WheelEvent("wheel", { bubbles: true, cancelable: true, deltaY: -100 }));
  });
  expect(zoom()).toBeCloseTo(1.331);
});

it("keeps a live drag when scrolling at the zoom limit leaves the scale unchanged", () => {
  const { svg, part, reorder } = setup();
  wheel(svg, -10000);
  fireEvent.mouseDown(part, { button: 0, clientX: 220, clientY: 80 });
  fireEvent.mouseMove(window, { clientX: 250, clientY: 80 });
  wheel(svg, -100);
  expect(screen.getByTestId("ring-drag-ghost")).toBeInTheDocument();
  fireEvent.mouseUp(part, { clientX: 250, clientY: 80 });
  expect(reorder).toHaveBeenCalledTimes(1);
});

it("removes its native wheel listener on unmount", () => {
  const { svg, unmount } = setup();
  expect(wheel(svg, -100).defaultPrevented).toBe(true);
  unmount();
  expect(wheel(svg, -100).defaultPrevented).toBe(false);
});

it("places two readable controls in the panel corner and removes zoom and rotation buttons", () => {
  const { container } = setup();
  expect(screen.queryByRole("button", { name: "缩小环图" })).not.toBeInTheDocument();
  expect(screen.queryByRole("button", { name: "放大环图" })).not.toBeInTheDocument();
  expect(screen.queryByRole("button", { name: "旋转环图" })).not.toBeInTheDocument();
  expect(screen.getByRole("button", { name: "重置环图大小" })).toHaveTextContent("重置大小");
  const toolbar = container.querySelector(".ring-controls")!;
  expect(toolbar.children).toHaveLength(2);
  expect(getComputedStyle(container.querySelector(".canvas")!).position).toBe("relative");
  expect(getComputedStyle(container.querySelector(".ring-wrap")!).position).toBe("static");
  expect(getComputedStyle(toolbar).position).toBe("absolute");
  for (const button of toolbar.querySelectorAll("button")) {
    expect(getComputedStyle(button).whiteSpace).toBe("nowrap");
    expect(getComputedStyle(button).width).toBe("auto");
  }
});
