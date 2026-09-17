import { afterEach, describe, expect, it, vi } from "vitest";
import { cleanup, fireEvent, render, screen } from "@testing-library/react";
import "@testing-library/jest-dom/vitest";
import "../src/style.css";
import { PlasmidRing } from "../src/PlasmidRing";
import type { ComponentType, Preview } from "../src/types";

const segment = (
  instance: string,
  id: string,
  type: ComponentType,
  start: number,
  end: number,
  kind: string = type === "t0" || type === "t1" ? "terminator" : type,
) => ({
  instance_id: instance,
  id,
  component_type: type,
  label: id,
  kind,
  start_bp: start,
  end_bp: end,
  length_bp: end - start + 1,
});

const preview: Preview = {
  valid: true,
  issues: [],
  warnings: [],
  sequence: "",
  length_bp: 1000,
  gc_percent: 50,
  sequence_sha256: "sequence",
  source_fingerprint: "source",
  manifest_revision: 1,
  segments: [
    segment("amp-a", "amp", "resistance", 101, 300),
    segment("amp-b", "amp", "resistance", 301, 500),
    segment("ori", "ori", "replication", 501, 700),
    segment("expression", "left_site", "expression", 701, 706, "restriction"),
    segment("expression", "expression", "expression", 707, 906),
    segment("expression", "right_site", "expression", 907, 912, "restriction"),
    segment("t1", "basic_seva_t1", "t1", 913, 942),
    segment("gap", "gap-14", "gap", 943, 956),
    segment("t0", "basic_seva_t0", "t0", 957, 981),
  ],
  features: [
    { label: "reverse gene", type: "CDS", kind: "gene", start_bp: 750, end_bp: 800, strand: -1 },
    { label: "EcoRI", type: "misc_feature", kind: "restriction", start_bp: 701, end_bp: 706, strand: 1 },
  ],
  enzymes: [],
};

function setup(data = preview, reorder = vi.fn()) {
  const props = {
    preview: data,
    construct: null,
    enzymes: [],
    componentOrder: ["amp-a", "amp-b", "ori", "expression", "t1", "gap", "t0"],
    onReorder: reorder,
  };
  const mounted = render(<PlasmidRing {...props} />);
  const svg = screen.getByRole("img", { name: "质粒环图" });
  vi.spyOn(svg, "getBoundingClientRect").mockReturnValue({
    left: 0, top: 0, width: 440, height: 440,
    right: 440, bottom: 440, x: 0, y: 0, toJSON() {},
  } as DOMRect);
  const part = (instance: string, kind?: string) => {
    const matches = Array.from(mounted.container.querySelectorAll(`[data-ring-instance="${instance}"]`));
    return matches.find((node) => !kind || node.querySelector("title")?.textContent?.startsWith(kind))!;
  };
  return { ...mounted, svg, part, props, reorder };
}

function click(part: Element, x = 220, y = 80) {
  fireEvent.mouseDown(part, { button: 0, clientX: x, clientY: y });
  fireEvent.mouseUp(part, { button: 0, clientX: x, clientY: y });
  fireEvent.click(part, { button: 0, detail: 1, clientX: x, clientY: y });
}

afterEach(() => { cleanup(); vi.restoreAllMocks(); });

describe("module reference-strand focus", () => {
  it("selects a small-movement click without reordering or changing DNA coordinates", () => {
    const { part, reorder } = setup();
    const amp = part("amp-a");
    const before = amp.querySelector("path")!.getAttribute("d");
    const otherBefore = part("amp-b").querySelector("path")!.getAttribute("d");
    fireEvent.mouseDown(amp, { button: 0, clientX: 220, clientY: 80 });
    fireEvent.mouseMove(window, { clientX: 223, clientY: 82 });
    fireEvent.mouseUp(amp, { clientX: 223, clientY: 82 });
    fireEvent.click(amp, { detail: 1, clientX: 223, clientY: 82 });
    expect(amp).toHaveAttribute("aria-pressed", "true");
    expect(screen.getByTestId("module-5-prime")).toHaveAttribute("data-bp", "101");
    expect(screen.getByTestId("module-3-prime")).toHaveAttribute("data-bp", "300");
    expect(amp.querySelector("path")).toHaveAttribute("data-start-bp", "101");
    expect(amp.querySelector("path")).toHaveAttribute("data-end-bp", "300");
    expect(amp.querySelector("path")!.getAttribute("d")).not.toBe(before);
    expect(part("amp-b").querySelector("path")).toHaveAttribute("d", otherBefore!);
    expect(screen.getByText("参考链方向 · 5′ → 3′")).toBeInTheDocument();
    expect(screen.queryByTestId("strand-arrow")).not.toBeInTheDocument();
    expect(reorder).not.toHaveBeenCalled();
  });

  it.each([
    ["amp-a", 101, 300], ["amp-b", 301, 500], ["ori", 501, 700],
    ["expression", 707, 906], ["t1", 913, 942], ["gap", 943, 956], ["t0", 957, 981],
  ])("marks the actual payload boundaries of %s", (instance, start, end) => {
    const { part } = setup();
    click(part(String(instance), instance === "expression" ? "expression" : undefined));
    expect(screen.getByTestId("module-5-prime")).toHaveAttribute("data-bp", String(start));
    expect(screen.getByTestId("module-3-prime")).toHaveAttribute("data-bp", String(end));
    expect(screen.getByTestId("module-direction")).toHaveAttribute("data-instance", String(instance));
  });

  it("switches independent duplicate instances and clears on repeated click, blank click or Escape", () => {
    const { part, svg } = setup();
    click(part("amp-a"));
    click(part("amp-b"));
    expect(part("amp-a")).toHaveAttribute("aria-pressed", "false");
    expect(part("amp-b")).toHaveAttribute("aria-pressed", "true");
    click(part("amp-b"));
    expect(screen.queryByTestId("module-direction")).not.toBeInTheDocument();
    click(part("amp-a"));
    fireEvent.click(svg);
    expect(screen.queryByTestId("module-direction")).not.toBeInTheDocument();
    click(part("amp-a"));
    fireEvent.keyDown(window, { key: "Escape" });
    expect(screen.queryByTestId("module-direction")).not.toBeInTheDocument();
  });

  it("keeps tiny gap labels separated and allows selection through its name leader", () => {
    const { container } = setup();
    click(container.querySelector(".gap-leader")!);
    const five = screen.getByTestId("module-5-prime-label");
    const three = screen.getByTestId("module-3-prime-label");
    expect(five).toHaveTextContent("5′");
    expect(three).toHaveTextContent("3′");
    expect(Math.abs(Number(five.getAttribute("y")) - Number(three.getAttribute("y")))).toBeGreaterThanOrEqual(24);
    expect(screen.getByTestId("module-5-prime").querySelector("path")).toHaveAttribute("fill", "none");
    expect(screen.getByTestId("module-3-prime").querySelector("path")).toHaveAttribute("fill", "none");
  });

  it("supports keyboard selection without mixing reference direction and reverse gene annotations", () => {
    const { part } = setup();
    const amp = part("amp-a");
    expect(amp).toHaveAttribute("role", "button");
    expect(amp).toHaveAttribute("tabindex", "0");
    fireEvent.keyDown(amp, { key: "Enter" });
    expect(screen.getByTestId("module-direction")).toBeInTheDocument();
    fireEvent.click(screen.getByRole("button", { name: "显示基因注释" }));
    expect(screen.getByTestId("strand-arrow")).toHaveAttribute("data-strand", "-1");
    expect(screen.getByTestId("module-5-prime")).toHaveAttribute("data-bp", "101");
    fireEvent.keyDown(amp, { key: " " });
    expect(screen.queryByTestId("module-direction")).not.toBeInTheDocument();
    expect(screen.getByTestId("strand-arrow")).toBeInTheDocument();
  });

  it("clears focus during actual dragging and suppresses the release click", () => {
    const { part, reorder } = setup();
    const amp = part("amp-a");
    click(amp);
    fireEvent.mouseDown(amp, { button: 0, clientX: 220, clientY: 80 });
    fireEvent.mouseMove(window, { clientX: 250, clientY: 80 });
    expect(screen.queryByTestId("module-direction")).not.toBeInTheDocument();
    expect(screen.getByTestId("ring-drag-ghost")).toBeInTheDocument();
    fireEvent.mouseUp(amp, { clientX: 250, clientY: 80 });
    fireEvent.click(amp, { detail: 1, clientX: 250, clientY: 80 });
    expect(reorder).toHaveBeenCalledTimes(1);
    expect(screen.queryByTestId("module-direction")).not.toBeInTheDocument();
    click(amp);
    expect(screen.getByTestId("module-direction")).toBeInTheDocument();
  });

  it("checks the movement threshold on release even without a mousemove event", () => {
    const { part, reorder } = setup();
    const amp = part("amp-a");
    fireEvent.mouseDown(amp, { button: 0, clientX: 220, clientY: 80 });
    fireEvent.mouseUp(amp, { clientX: 225, clientY: 80 });
    fireEvent.click(amp, { detail: 1, clientX: 225, clientY: 80 });
    expect(reorder).toHaveBeenCalledTimes(1);
    expect(screen.queryByTestId("module-direction")).not.toBeInTheDocument();
  });

  it.each(["Escape", "blur", "geometry"])("does not select after a gesture cancelled by %s", (reason) => {
    const { part, props, rerender, reorder } = setup();
    const amp = part("amp-a");
    fireEvent.mouseDown(amp, { button: 0, clientX: 220, clientY: 80 });
    if (reason === "Escape") fireEvent.keyDown(window, { key: "Escape" });
    else if (reason === "blur") fireEvent.blur(window);
    else rerender(<PlasmidRing {...props} preview={{ ...preview, length_bp: 1001 }} />);
    fireEvent.mouseUp(amp, { clientX: 220, clientY: 80 });
    fireEvent.click(amp, { detail: 1, clientX: 220, clientY: 80 });
    expect(reorder).not.toHaveBeenCalled();
    expect(screen.queryByTestId("module-direction")).not.toBeInTheDocument();
    click(amp);
    expect(screen.getByTestId("module-direction")).toBeInTheDocument();
  });

  it("preserves selection through identical snapshots, coordinate updates, wheel zoom and reset", () => {
    const { part, props, rerender } = setup();
    click(part("amp-a"));
    rerender(<PlasmidRing {...props} preview={{ ...preview, segments: preview.segments.map((p) => ({ ...p })) }} />);
    expect(screen.getByTestId("module-direction")).toBeInTheDocument();
    const moved = { ...preview, segments: preview.segments.map((p) => p.instance_id === "amp-a" ? { ...p, start_bp: 201, end_bp: 400 } : p) };
    rerender(<PlasmidRing {...props} preview={moved} />);
    expect(screen.getByTestId("module-5-prime")).toHaveAttribute("data-bp", "201");
    expect(screen.getByTestId("module-3-prime")).toHaveAttribute("data-bp", "400");
    fireEvent.wheel(screen.getByRole("img", { name: "质粒环图" }), { deltaY: -100 });
    const transform = screen.getByTestId("module-direction").parentElement!.getAttribute("transform");
    expect(transform).toContain("scale(1.1)");
    fireEvent.click(screen.getByRole("button", { name: "重置环图大小" }));
    expect(screen.getByTestId("module-direction").parentElement!.getAttribute("transform")).toContain("scale(1)");
    expect(screen.getByTestId("module-5-prime")).toHaveAttribute("data-bp", "201");
  });

  it("keeps expression selected when provisional and validated payload IDs differ", () => {
    const provisional = {
      ...preview,
      source_fingerprint: "",
      segments: preview.segments
        .filter((p) => p.kind !== "restriction")
        .map((p) => p.kind === "expression" ? { ...p, id: "1", start_bp: 701, end_bp: 900 } : p),
    };
    const { part, props, rerender } = setup(provisional);
    click(part("expression"));
    expect(screen.getByTestId("module-5-prime")).toHaveAttribute("data-bp", "701");
    rerender(<PlasmidRing {...props} preview={preview} />);
    expect(screen.getByTestId("module-5-prime")).toHaveAttribute("data-bp", "707");
    expect(screen.getByTestId("module-3-prime")).toHaveAttribute("data-bp", "906");
    rerender(<PlasmidRing {...props} preview={provisional} />);
    expect(screen.getByTestId("module-5-prime")).toHaveAttribute("data-bp", "701");
    expect(screen.getByTestId("module-3-prime")).toHaveAttribute("data-bp", "900");
  });

  it.each(["blur", "geometry"])("accepts a fresh blank click after %s without a release click", (reason) => {
    const { part, props, rerender, svg } = setup();
    click(part("amp-a"));
    fireEvent.mouseDown(part("amp-a"), { button: 0, clientX: 220, clientY: 80 });
    if (reason === "blur") fireEvent.blur(window);
    else rerender(<PlasmidRing {...props} preview={{ ...preview, length_bp: 1001 }} />);
    fireEvent.mouseUp(window, { clientX: 220, clientY: 80 });
    expect(screen.getByTestId("module-direction")).toBeInTheDocument();
    click(svg, 220, 220);
    expect(screen.queryByTestId("module-direction")).not.toBeInTheDocument();
  });

  it.each(["delete", "replace", "source"])("clears stale selection when the module changes: %s", (reason) => {
    const { part, props, rerender } = setup();
    click(part("amp-a"));
    expect(screen.getByTestId("module-direction")).toBeInTheDocument();
    const changed = reason === "source"
      ? { ...preview, source_fingerprint: "new-source" }
      : { ...preview, segments: preview.segments.flatMap((p) => p.instance_id !== "amp-a" ? [p] : reason === "delete" ? [] : [{ ...p, id: "kan" }]) };
    rerender(<PlasmidRing {...props} preview={changed} />);
    expect(screen.queryByTestId("module-direction")).not.toBeInTheDocument();
    rerender(<PlasmidRing {...props} />);
    expect(screen.queryByTestId("module-direction")).not.toBeInTheDocument();
  });

  it("can focus a single full-circle construct without a reorder callback", () => {
    const { container } = render(<PlasmidRing preview={null} construct={{
      id: 1, name: "Source", length_bp: 1000, gc_percent: 50, sequence_sha256: "source", features: [],
    }} enzymes={[]} />);
    const circle = screen.getByTestId("full-segment");
    const before = Number(circle.getAttribute("r"));
    click(container.querySelector('[data-ring-component="expression"]')!);
    expect(Number(circle.getAttribute("r"))).toBeGreaterThan(before);
    expect(screen.getByTestId("module-5-prime")).toHaveAttribute("data-bp", "1");
    expect(screen.getByTestId("module-3-prime")).toHaveAttribute("data-bp", "1000");
    for (const end of ["module-5-prime", "module-3-prime"]) {
      const dot = screen.getByTestId(end).querySelector("circle")!;
      expect(Number(dot.getAttribute("cx"))).toBeCloseTo(220);
      expect(Number(dot.getAttribute("cy"))).toBeLessThan(80);
    }
    const five = screen.getByTestId("module-5-prime-label");
    const three = screen.getByTestId("module-3-prime-label");
    expect(Math.hypot(
      Number(five.getAttribute("x")) - Number(three.getAttribute("x")),
      Number(five.getAttribute("y")) - Number(three.getAttribute("y")),
    )).toBeGreaterThanOrEqual(24);
  });
});
