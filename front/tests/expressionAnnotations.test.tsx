import { afterEach, expect, it, vi } from "vitest";
import { cleanup, fireEvent, render, screen, within } from "@testing-library/react";
import "@testing-library/jest-dom/vitest";
import "../src/style.css";
import "../src/controls.css";
import { PlasmidRing } from "../src/PlasmidRing";
import { reorderGeometry } from "../src/componentOrder";
import type { Feature, Preview } from "../src/types";

const feature = (label: string, type: string, start: number, end: number, strand = 1): Feature => ({
  label, type, start_bp: start, end_bp: end, strand,
  kind: type === "CDS" ? "gene" : "expression",
});
const features = [
  feature("expression_cassette_1", "misc_feature", 1, 2525),
  feature("expression_cassette_2", "misc_feature", 2526, 3553),
  feature("cassette_1_promoter", "promoter", 1, 35),
  feature("cassette_1_P21683_RBS", "RBS", 36, 47),
  feature("cassette_1_P21683_CDS", "CDS", 48, 977),
  feature("cassette_1_P21685_RBS", "RBS", 978, 989),
  feature("cassette_1_P21685_CDS", "CDS", 990, 2468),
  feature("cassette_1_terminator", "terminator", 2469, 2525),
  feature("cassette_2_promoter", "promoter", 2526, 2560),
  feature("cassette_2_P22873_RBS", "RBS", 2561, 2572),
  feature("cassette_2_P22873_CDS", "CDS", 2573, 3496),
  feature("cassette_2_terminator", "terminator", 3497, 3553),
];
const preview: Preview = {
  valid: true, issues: [], warnings: [], sequence: "", length_bp: 6138,
  gc_percent: 50, sequence_sha256: "sequence", source_fingerprint: "source",
  manifest_revision: 1, enzymes: [], features,
  segments: [
    { id: "expression", label: "expression", kind: "expression", component_type: "expression", instance_id: "expression", start_bp: 1, end_bp: 3553, length_bp: 3553 },
    { id: "ori", label: "ori", kind: "replication", component_type: "replication", instance_id: "ori", start_bp: 3554, end_bp: 6138, length_bp: 2585 },
  ],
};

function setup(data = preview) {
  const reorder = vi.fn();
  const props = { preview: data, construct: null, enzymes: [], componentOrder: ["expression", "ori"], onReorder: reorder };
  const mounted = render(<section className="canvas"><PlasmidRing {...props} /></section>);
  const svg = screen.getByRole("img", { name: "质粒环图" });
  const show = () => fireEvent.click(screen.getByRole("button", { name: "显示基因注释" }));
  const nodes = (role: string) => svg.querySelectorAll(`[data-annotation-role="${role}"]`);
  return { ...mounted, props, svg, show, nodes, reorder };
}

afterEach(() => { cleanup(); vi.restoreAllMocks(); });

it("draws all ten internal elements and two cassette groups only with annotations enabled", () => {
  const { svg, show, nodes } = setup();
  expect(nodes("gene")).toHaveLength(0);
  show();
  expect(nodes("promoter")).toHaveLength(2);
  expect(nodes("rbs")).toHaveLength(3);
  expect(nodes("gene")).toHaveLength(3);
  expect(nodes("terminator")).toHaveLength(2);
  expect(svg.querySelectorAll("[data-cassette-group]")).toHaveLength(2);
  expect(svg.querySelectorAll("[data-gene-label]")).toHaveLength(3);
  expect(Array.from(svg.querySelectorAll("[data-gene-label] textPath")).map(n => n.textContent)).toEqual(["P21683", "P21685", "P22873"]);
  expect(screen.getByTestId("annotation-key")).toHaveTextContent("启动子");
  expect(screen.getByTestId("annotation-key")).toHaveTextContent("RBS");
  fireEvent.click(screen.getByRole("button", { name: "隐藏基因注释" }));
  expect(nodes("gene")).toHaveLength(0);
});

it("uses feature types for colors and separates adjacent promoter and tiny RBS markers", () => {
  const { nodes, show } = setup();
  show();
  const promoter = nodes("promoter")[0];
  const rbs = nodes("rbs")[0];
  expect(promoter).toHaveAttribute("data-start-bp", "1");
  expect(promoter).toHaveAttribute("data-end-bp", "35");
  expect(rbs).toHaveAttribute("data-start-bp", "36");
  expect(rbs).toHaveAttribute("data-end-bp", "47");
  expect(promoter.getAttribute("data-track-radius")).not.toBe(rbs.getAttribute("data-track-radius"));
  expect(promoter.querySelector(".annotation-symbol")).toHaveTextContent("P");
  expect(rbs.querySelector(".annotation-symbol")).toHaveTextContent("R");
  expect(nodes("terminator")[0].querySelector(".annotation-symbol")).toHaveTextContent("T");
  expect(getComputedStyle(promoter.querySelector(".annotation-symbol")!).fill)
    .not.toBe(getComputedStyle(rbs.querySelector(".annotation-symbol")!).fill);
});

it("keeps full labels, lengths and strand in hover text and organizes details by cassette", () => {
  const { nodes, container, show } = setup();
  show();
  expect(nodes("gene")[0].querySelector("title")).toHaveTextContent("cassette_1_P21683_CDS");
  expect(nodes("gene")[0].querySelector("title")).toHaveTextContent("930 bp");
  expect(nodes("gene")[0].querySelector("title")).toHaveTextContent("48–977");
  const groups = container.querySelectorAll("[data-cassette-details]");
  expect(groups).toHaveLength(2);
  expect(within(groups[0] as HTMLElement).getByText("表达盒1")).toBeInTheDocument();
  expect(groups[0]).toHaveTextContent("P21683");
  expect(groups[0]).toHaveTextContent("P21685");
  expect(groups[1]).toHaveTextContent("P22873");
});

it.each([1, -1, 0])("connects a CDS arrow inside its actual boundaries for strand %s", (strand) => {
  const { svg, show, nodes } = setup({ ...preview, length_bp: 1000, features: [feature("gene", "CDS", 251, 500, strand)], segments: [{ ...preview.segments[0], end_bp: 1000, length_bp: 1000 }] });
  show();
  const path = nodes("gene")[0].querySelector("path[data-start-bp]")!;
  expect(path).toHaveAttribute("data-start-bp", "251");
  expect(path).toHaveAttribute("data-end-bp", "500");
  expect(svg.querySelector("polygon[data-testid='strand-arrow']")).toBeNull();
  if (strand === 1) expect(path.getAttribute("d")).toContain("L 220 323");
  if (strand === -1) expect(path.getAttribute("d")).toMatch(/^M 323 220/);
  if (strand === 0) {
    expect(svg.querySelector("[data-testid='strand-arrow']")).toBeNull();
    expect(nodes("gene")[0].querySelector("title")).toHaveTextContent("方向未标注");
  } else expect(path).toHaveAttribute("data-strand", String(strand));
});

it("displays short and long-name genes with safe leaders and retains unabridged hover information", () => {
  const longName = "cassette_1_VERY_LONG_PROTEIN_IDENTIFIER_THAT_CANNOT_FIT_CDS";
  const { svg, show, nodes } = setup({ ...preview, features: [feature(longName, "CDS", 200, 213)] });
  show();
  expect(nodes("gene")).toHaveLength(1);
  const label = svg.querySelector("[data-gene-label]")!;
  expect(label).toBeInTheDocument();
  expect(label).toHaveAttribute("data-label-placement", "leader");
  expect(label.querySelector("title")).toHaveTextContent(longName);
  expect(getComputedStyle(label).pointerEvents).toBe("auto");
  expect(nodes("gene")[0].querySelector("title")).toHaveTextContent(longName);
  expect(nodes("gene")[0].querySelector("path")!.getAttribute("d")).not.toMatch(/NaN|Infinity/);
});

it("does not infer a gene from generic wrappers with a misleading kind", () => {
  const { nodes, show } = setup({ ...preview, features: [
    ...features, { ...feature("generic wrapper", "misc_feature", 1, 3553), kind: "gene" },
  ] });
  show();
  expect(nodes("gene")).toHaveLength(3);
  const wrapper = screen.getByText("generic wrapper", { selector: "li > span" }).closest("li")!;
  const gene = screen.getByText("P21683", { selector: "li > span" }).closest("li")!;
  expect(getComputedStyle(wrapper.querySelector("i")!).backgroundColor)
    .not.toBe(getComputedStyle(gene.querySelector("i")!).backgroundColor);
});

it("keeps non-expression CDS but excludes outside terminators and module wrappers from internal elements", () => {
  const { nodes, svg, show } = setup({ ...preview, features: [
    ...features,
    feature("Rep101", "CDS", 3700, 4650),
    feature("T1", "terminator", 5700, 5804),
    feature("full module", "misc_feature", 1, 3553),
  ] });
  show();
  expect(nodes("gene")).toHaveLength(4);
  expect(nodes("terminator")).toHaveLength(2);
  expect(svg.querySelectorAll("[data-cassette-group]")).toHaveLength(2);
});

it("shows flat details when no cassette annotation exists and handles a full-circle CDS", () => {
  const { svg, show, container, nodes } = setup({ ...preview, length_bp: 1000, features: [feature("whole gene", "CDS", 1, 1000, 0)], segments: [{ ...preview.segments[0], end_bp: 1000, length_bp: 1000 }] });
  show();
  expect(container.querySelectorAll("[data-cassette-details]")).toHaveLength(0);
  expect(svg.querySelectorAll("[data-cassette-group]")).toHaveLength(0);
  expect(nodes("gene")[0].querySelector("path")!.getAttribute("d")).not.toMatch(/NaN|Infinity/);
  expect((nodes("gene")[0].querySelector("path")!.getAttribute("d")!.match(/A /g) || []).length).toBeGreaterThanOrEqual(4);
});

it("carries internal positions through reordering and preserves focus and annotation state on wheel/reset", () => {
  const { props, rerender, svg, show, nodes, reorder } = setup();
  show();
  const moved = reorderGeometry(preview, ["ori", "expression"]);
  rerender(<section className="canvas"><PlasmidRing {...props} preview={moved} componentOrder={["ori", "expression"]} /></section>);
  expect(nodes("promoter")[0]).toHaveAttribute("data-start-bp", "2586");
  expect(nodes("gene")[0]).toHaveAttribute("data-start-bp", "2633");
  expect(nodes("gene")[0]).toHaveAttribute("data-end-bp", "3562");
  fireEvent.click(svg.querySelector('[data-ring-instance="expression"]')!);
  const direction = screen.getByTestId("module-direction");
  fireEvent.wheel(svg, { deltaY: -10000 });
  expect(nodes("gene")[0].closest("g[transform]")).toBe(direction.parentElement);
  expect(direction.parentElement!.getAttribute("transform")).toContain("scale(2)");
  fireEvent.click(screen.getByRole("button", { name: "重置环图大小" }));
  expect(screen.getByTestId("module-direction")).toBeInTheDocument();
  expect(nodes("gene")).toHaveLength(3);
  expect(reorder).not.toHaveBeenCalled();
});

it("keeps labels for distinct SVG instances bound to their own paths", () => {
  const props = { preview, construct: null, enzymes: [] };
  const { container } = render(<><PlasmidRing {...props} /><PlasmidRing {...props} /></>);
  screen.getAllByRole("button", { name: "显示基因注释" }).forEach(button => fireEvent.click(button));
  const ids = Array.from(container.querySelectorAll("defs path[id]")).map(n => n.id);
  expect(ids.length).toBeGreaterThanOrEqual(6);
  expect(new Set(ids).size).toBe(ids.length);
  for (const svg of container.querySelectorAll("svg")) {
    for (const text of svg.querySelectorAll("textPath")) {
      const id = text.getAttribute("href")!.slice(1);
      expect(Array.from(svg.querySelectorAll("[id]")).some(n => n.id === id)).toBe(true);
    }
  }
});
