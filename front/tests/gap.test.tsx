import { afterEach, beforeEach, expect, it, vi } from "vitest";
import { cleanup, fireEvent, render, screen, waitFor, within } from "@testing-library/react";
import "@testing-library/jest-dom/vitest";
import { Workbench } from "../src/Workbench";
import { restoreComponents } from "../src/componentSelection";
import { PlasmidRing } from "../src/PlasmidRing";
import type { Module } from "../src/types";
import "../src/style.css";

const module = (id: string, type: Module["type"], length = 4, extra = {}) => ({
  id,
  name: id,
  type,
  sequence: "ATGC".repeat(Math.ceil(length / 4)).slice(0, length),
  length_bp: length,
  gc_percent: 50,
  notes: [],
  ...extra,
});
const library = {
  resistance: [module("amp", "resistance")],
  replication: [module("ori", "replication")],
  terminator: [],
  gap: [
    module("gap-r", "gap", 14, { aliases: ["resistance_to_replication"], purpose: "native interval", evidence_status: "source_sequence_verified" }),
    module("gap-o", "gap", 14, { aliases: ["replication_to_t1"], purpose: "native interval", evidence_status: "verified" }),
    module("landing", "gap", 14, { aliases: ["landing_pad_spacer"], purpose: "neutral spacer", evidence_status: "candidate" }),
  ],
};
const context = {
  project: { target: "gap-demo", name: "Gap demo", host: "E. coli", manifest_revision: 1, source_fingerprint: "source" },
  ready: true,
  issues: [],
  construct: { id: 1, name: "Source", length_bp: 100, gc_percent: 50, sequence_sha256: "source", features: [] },
  restriction_enzymes: [],
  selection: null,
  result: null,
  active_job: null,
};

function install() {
  const fetcher = vi.fn((url: string, init?: RequestInit) => {
    const body = JSON.parse(String(init?.body || "{}"));
    let cursor = 1;
    const segments = (body.components || []).map((component: { instance_id: string; component_type: string; module_id?: string }) => {
      const length = component.component_type === "expression"
        ? 100
        : Object.values(library)
            .flat()
            .find((item) => item.id === component.module_id)?.length_bp || 0;
      const segment = {
        id: component.module_id || "expression",
        instance_id: component.instance_id,
        component_type: component.component_type,
        label: component.module_id || "Source",
        kind: component.component_type,
        start_bp: cursor,
        end_bp: cursor + length - 1,
        length_bp: length,
      };
      cursor += length;
      return segment;
    });
    return Promise.resolve(new Response(JSON.stringify(
      url === "/api/modules" ? library : url === "/api/preview" ? {
        ...body, valid: true, issues: [], warnings: [], terminator_warnings: [], sequence: "", length_bp: cursor - 1,
        gc_percent: 50, sequence_sha256: "preview", source_fingerprint: "source", manifest_revision: 1, segments, features: [], enzymes: [],
      } : context,
    )));
  });
  vi.stubGlobal("fetch", fetcher);
  return fetcher;
}

beforeEach(() => { localStorage.clear(); install(); });
afterEach(() => { cleanup(); vi.restoreAllMocks(); vi.unstubAllGlobals(); });

it("adds three independent gaps, removes only the middle instance, and previews v2", async () => {
  const fetcher = install();
  const { container } = render(<Workbench />);
  fireEvent.click(await screen.findByRole("button", { name: "选择 amp" }));
  fireEvent.click(screen.getByRole("button", { name: "选择 ori" }));
  const gap = await screen.findByRole("button", { name: "选择 landing" });
  fireEvent.click(gap);
  fireEvent.click(gap);
  fireEvent.click(gap);
  const rows = await screen.findAllByTestId("gap-slot");
  const ids = rows.map((row) => row.closest(".component-row")!.getAttribute("data-component-instance"));
  expect(new Set(ids).size).toBe(3);
  expect(container.querySelectorAll('[data-ring-component="gap"]')).toHaveLength(3);
  fireEvent.click(within(rows[1]).getByRole("button", { name: "清空间隔区段" }));
  expect(screen.getAllByTestId("gap-slot").map((row) => row.closest(".component-row")!.getAttribute("data-component-instance"))).toEqual([ids[0], ids[2]]);
  await waitFor(() => expect(fetcher.mock.calls.some(([url, init]) =>
    url === "/api/preview" && JSON.parse(String(init?.body)).assembly_schema_version === 2,
  )).toBe(true));
});

it("keeps v2 gapless selections gapless and expands legacy owners with deterministic collision-safe IDs", () => {
  const v2 = restoreComponents({ assembly_schema_version: 2, components: [
    { instance_id: "replication", component_type: "replication", module_id: "ori" },
    { instance_id: "expression", component_type: "expression" },
    { instance_id: "resistance", component_type: "resistance", module_id: "amp" },
  ] }, library);
  expect(v2.map((component) => component.component_type)).not.toContain("gap");

  const legacy = restoreComponents({ components: [
    { instance_id: "legacy-gap-0-resistance", component_type: "resistance", module_id: "amp" },
    { instance_id: "replication", component_type: "replication", module_id: "ori" },
    { instance_id: "expression", component_type: "expression" },
  ] }, library);
  expect(legacy).toEqual([
    { instance_id: "legacy-gap-0-resistance", component_type: "resistance", module_id: "amp" },
    { instance_id: "legacy-gap-0-resistance-1", component_type: "gap", module_id: "gap-r" },
    { instance_id: "replication", component_type: "replication", module_id: "ori" },
    { instance_id: "legacy-gap-1-replication", component_type: "gap", module_id: "gap-o" },
    { instance_id: "expression", component_type: "expression" },
  ]);
});

it("keeps an explicit unversioned gap and gives a 14 bp gap a grey coordinate callout", () => {
  const restored = restoreComponents({ components: [
    { instance_id: "replication", component_type: "replication", module_id: "ori" },
    { instance_id: "gap-instance", component_type: "gap", module_id: "landing" },
    { instance_id: "expression", component_type: "expression" },
    { instance_id: "resistance", component_type: "resistance", module_id: "amp" },
  ] }, library);
  expect(restored.filter((component) => component.component_type === "gap")).toEqual([
    { instance_id: "gap-instance", component_type: "gap", module_id: "landing" },
  ]);

  const { container } = render(
    <PlasmidRing
      preview={{
        valid: true, issues: [], warnings: [], sequence: "", length_bp: 1000,
        gc_percent: 50, sequence_sha256: "preview", source_fingerprint: "source", manifest_revision: 1,
        components: restored,
        segments: [{ id: "landing", instance_id: "gap-instance", component_type: "gap", label: "Landing spacer", kind: "gap", start_bp: 100, end_bp: 113, length_bp: 14 }],
        features: [], enzymes: [],
      }}
      construct={null}
      enzymes={[]}
    />,
  );
  expect(container.querySelector(".gap-callout")).toHaveTextContent("Landing spacer");
  expect(container.querySelector('path[data-start-bp="100"]')).toHaveAttribute("data-end-bp", "113");
  expect(container.querySelector(".gap-leader path")).toHaveAttribute("stroke", "#8b949e");
});

it("uses selected scalar positions, reserves filtered legacy IDs, and skips an empty replication placeholder", () => {
  const scalar = restoreComponents({
    component_order: ["t1", "expression", "t0", "replication", "resistance"],
    t1_id: null,
    t0_id: null,
    replication_id: "ori",
    resistance_id: "amp",
  }, library);
  expect(scalar.map((component) => component.instance_id)).toEqual([
    "expression",
    "replication",
    "legacy-gap-1-replication",
    "resistance",
    "legacy-gap-2-resistance",
  ]);

  const collision = restoreComponents({ components: [
    { instance_id: "legacy-gap-1-replication", component_type: "t0", module_id: "removed" },
    { instance_id: "replication", component_type: "replication", module_id: "ori" },
    { instance_id: "expression", component_type: "expression" },
    { instance_id: "resistance", component_type: "resistance", module_id: "amp" },
  ] }, library);
  expect(collision.map((component) => component.instance_id)).toContain("legacy-gap-1-replication-1");

  const placeholder = restoreComponents({ components: [
    { instance_id: "replication", component_type: "replication" },
    { instance_id: "expression", component_type: "expression" },
    { instance_id: "resistance", component_type: "resistance", module_id: "amp" },
  ] }, library);
  expect(placeholder).toEqual([
    { instance_id: "replication", component_type: "replication" },
    { instance_id: "expression", component_type: "expression" },
    { instance_id: "resistance", component_type: "resistance", module_id: "amp" },
    { instance_id: "legacy-gap-2-resistance", component_type: "gap", module_id: "gap-r" },
  ]);
});

it("shows a human-readable evidence label for a selected gap", async () => {
  render(<Workbench />);
  fireEvent.click(await screen.findByRole("button", { name: "选择 gap-r" }));
  fireEvent.click(screen.getByText("组件详情与 DNA"));
  expect(await screen.findByText(/来源序列已核对/)).toBeInTheDocument();
  expect(screen.queryByText("source_sequence_verified")).not.toBeInTheDocument();
});

it("applies the grey gap color to computed callout text and the drag insertion marker", () => {
  const { container } = render(
    <PlasmidRing
      preview={{
        valid: true, issues: [], warnings: [], sequence: "", length_bp: 114,
        gc_percent: 50, sequence_sha256: "preview", source_fingerprint: "source", manifest_revision: 1,
        components: [
          { instance_id: "gap", component_type: "gap", module_id: "landing" },
          { instance_id: "expression", component_type: "expression" },
        ],
        segments: [
          { id: "landing", instance_id: "gap", component_type: "gap", label: "Landing spacer", kind: "gap", start_bp: 1, end_bp: 14, length_bp: 14 },
          { id: "source", instance_id: "expression", component_type: "expression", label: "Source", kind: "expression", start_bp: 15, end_bp: 114, length_bp: 100 },
        ],
        features: [], enzymes: [],
      }}
      construct={null}
      enzymes={[]}
      componentOrder={["gap", "expression"]}
      onReorder={() => {}}
    />,
  );
  expect(getComputedStyle(container.querySelector(".gap-callout")!).fill).toBe("rgb(139, 148, 158)");
  const svg = screen.getByRole("img", { name: "质粒环图" });
  vi.spyOn(svg, "getBoundingClientRect").mockReturnValue({
    x: 0, y: 0, top: 0, left: 0, bottom: 440, right: 440, width: 440, height: 440,
    toJSON: () => ({}),
  });
  fireEvent.mouseDown(container.querySelector('[data-ring-instance="gap"]')!, { button: 0, clientX: 220, clientY: 68 });
  fireEvent.mouseMove(window, { clientX: 370, clientY: 220 });
  const marker = screen.getByTestId("ring-insertion-marker");
  expect(getComputedStyle(marker).stroke).toBe("rgb(139, 148, 158)");
  expect(getComputedStyle(marker).strokeWidth).toBe("4");
});

it("persists a reordered duplicate gap list across reload and hides a stale ordered result", async () => {
  const saved = [
    { instance_id: "replication", component_type: "replication", module_id: "ori" },
    { instance_id: "expression", component_type: "expression" },
    { instance_id: "gap-a", component_type: "gap", module_id: "landing" },
    { instance_id: "gap-b", component_type: "gap", module_id: "landing" },
    { instance_id: "resistance", component_type: "resistance", module_id: "amp" },
  ];
  localStorage.setItem("plasmid:gap-demo", JSON.stringify({ assembly_schema_version: 2, components: saved }));
  const staleResult = {
    id: "stale", components: saved, source_fingerprint: "source", manifest_revision: 1,
    files: [{ id: "final_genbank", label: "GenBank", filename: "stale.gb", format: "gb", url: "/stale.gb" }], warnings: [],
  };
  const fetcher = vi.fn((url: string, init?: RequestInit) => {
    const body = JSON.parse(String(init?.body || "{}"));
    if (url === "/api/modules") return Promise.resolve(new Response(JSON.stringify(library)));
    if (url === "/api/preview") return Promise.resolve(new Response(JSON.stringify({
      ...body, valid: true, issues: [], warnings: [], terminator_warnings: [], sequence: "", length_bp: 136,
      gc_percent: 50, sequence_sha256: "preview", source_fingerprint: "source", manifest_revision: 1,
      segments: body.components.map((component: { instance_id: string; component_type: string; module_id?: string }, index: number) => ({
        id: component.module_id || "source", instance_id: component.instance_id, component_type: component.component_type,
        label: component.module_id || "Source", kind: component.component_type, start_bp: index * 14 + 1, end_bp: index * 14 + 14, length_bp: 14,
      })), features: [], enzymes: [],
    })));
    return Promise.resolve(new Response(JSON.stringify({ ...context, result: staleResult })));
  });
  vi.stubGlobal("fetch", fetcher);
  const { container, unmount } = render(<Workbench />);
  await screen.findByRole("link", { name: "下载 GenBank" });
  const rows = Array.from(container.querySelectorAll<HTMLElement>(".component-row"));
  rows.forEach((row, index) => vi.spyOn(row, "getBoundingClientRect").mockReturnValue({
    x: 0, y: index * 60, top: index * 60, left: 0, right: 200, bottom: index * 60 + 50, width: 200, height: 50, toJSON: () => ({}),
  }));
  fireEvent.mouseDown(screen.getAllByTestId("gap-slot")[1], { button: 0, clientX: 20, clientY: 190 });
  fireEvent.mouseMove(window, { clientX: 20, clientY: 0 });
  fireEvent.mouseUp(window, { clientX: 20, clientY: 0 });
  const expected = ["gap-b", "replication", "expression", "gap-a", "resistance"];
  await waitFor(() => expect(JSON.parse(localStorage.getItem("plasmid:gap-demo")!).components.map((component: { instance_id: string }) => component.instance_id)).toEqual(expected));
  expect(screen.queryByRole("link", { name: "下载 GenBank" })).not.toBeInTheDocument();
  unmount();
  render(<Workbench />);
  await waitFor(() => expect(Array.from(document.querySelectorAll("[data-component-instance]")).map((row) => row.getAttribute("data-component-instance"))).toEqual(expected));
  expect(screen.queryByRole("link", { name: "下载 GenBank" })).not.toBeInTheDocument();
});
