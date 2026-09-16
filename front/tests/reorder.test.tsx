import { afterEach, beforeEach, describe, expect, it, vi } from "vitest";
import {
  act,
  cleanup,
  fireEvent,
  render,
  screen,
  waitFor,
} from "@testing-library/react";
import "@testing-library/jest-dom/vitest";
import { Workbench } from "../src/Workbench";
import { PlasmidRing } from "../src/PlasmidRing";
import type { Preview, ComponentOrder } from "../src/types";

const defaultOrder: ComponentOrder = [
  "resistance",
  "replication",
  "t1",
  "expression",
  "t0",
];
const movedOrder: ComponentOrder = [
  "replication",
  "t1",
  "expression",
  "t0",
  "resistance",
];
const modules = {
  resistance: [
    {
      id: "amp",
      name: "AmpR",
      type: "resistance",
      sequence: "ATGC",
      length_bp: 20,
      gc_percent: 50,
      notes: [],
    },
    {
      id: "kan",
      name: "KanR",
      type: "resistance",
      sequence: "ATGC",
      length_bp: 20,
      gc_percent: 50,
      notes: [],
    },
  ],
  replication: [
    {
      id: "ori",
      name: "pUC ori",
      type: "replication",
      sequence: "ATGC",
      length_bp: 30,
      gc_percent: 50,
      notes: [],
    },
  ],
};
const context = {
  project: {
    target: "drag-demo",
    name: "Demo",
    host: "E. coli",
    manifest_revision: 1,
    source_fingerprint: "source",
  },
  ready: true,
  issues: [],
  construct: {
    id: 1,
    name: "Source",
    length_bp: 40,
    gc_percent: 50,
    sequence_sha256: "source",
    features: [],
  },
  restriction_enzymes: [],
  selection: {
    resistance_id: "amp",
    replication_id: "ori",
    component_order: defaultOrder,
  },
  active_job: null,
  result: null,
};
let cursor = 1;
const segments = [
  ["amp", "resistance", 20],
  ["module_interval", "linker", 2],
  ["ori", "replication", 30],
  ["t1_interval", "linker", 3],
  ["t1", "terminator", 5],
  ["left_site", "restriction", 6],
  ["source", "expression", 40],
  ["right_site", "restriction", 6],
].map(([id, kind, size]) => {
  const start = cursor;
  cursor += Number(size);
  return {
    id: String(id),
    label: String(id),
    kind: String(kind),
    length_bp: Number(size),
    start_bp: start,
    end_bp: cursor - 1,
  };
});
const preview: Preview = {
  valid: true,
  issues: [],
  warnings: [],
  sequence: "",
  length_bp: 112,
  gc_percent: 50,
  sequence_sha256: "preview",
  source_fingerprint: "source",
  manifest_revision: 1,
  segments,
  features: [
    {
      label: "EcoRI",
      type: "misc_feature",
      kind: "restriction",
      strand: 1,
      start_bp: 107,
      end_bp: 112,
    },
  ],
  enzymes: [],
};
const result = {
  id: "old",
  resistance_id: "amp",
  replication_id: "ori",
  component_order: defaultOrder,
  source_fingerprint: "source",
  manifest_revision: 1,
  files: [
    {
      id: "final_genbank",
      label: "GenBank",
      filename: "old.gb",
      url: "/old.gb",
    },
  ],
  warnings: [],
};
const response = (body: unknown) =>
  Promise.resolve(new Response(JSON.stringify(body)));
function install(
  readPreview: (init?: RequestInit) => Promise<Response> = () =>
    response({ ...preview, component_order: defaultOrder }),
  source: unknown = context,
) {
  const mock = vi.fn((url: string, init?: RequestInit) =>
    url.includes("modules")
      ? response(modules)
      : url.includes("preview")
        ? readPreview(init)
        : url.includes("generate")
          ? response({
              id: "j",
              status: "succeeded",
              stage: "done",
              message: "done",
              result: {
                ...result,
                component_order: JSON.parse(String(init?.body)).component_order,
              },
            })
          : response(source),
  );
  vi.stubGlobal("fetch", mock);
  return mock;
}
function rect(left: number, top: number, width: number, height: number) {
  return {
    left,
    top,
    width,
    height,
    right: left + width,
    bottom: top + height,
    x: left,
    y: top,
    toJSON() {},
  } as DOMRect;
}
function layoutRows() {
  ["resistance", "replication", "expression"].forEach((type, index) => {
    const row = document.querySelector(`[data-component-type='${type}']`)!;
    expect(row).toBeInTheDocument();
    vi.spyOn(row, "getBoundingClientRect").mockReturnValue(
      rect(600, 100 + index * 60, 200, 50),
    );
  });
}
function listDrag() {
  layoutRows();
  fireEvent.mouseDown(screen.getByTestId("resistance-slot"), {
    button: 0,
    clientX: 620,
    clientY: 120,
  });
  fireEvent.mouseMove(window, { clientX: 620, clientY: 290 });
  fireEvent.mouseUp(window, { clientX: 620, clientY: 290 });
}
function rowOrder() {
  return Array.from(document.querySelectorAll("[data-component-type]")).map(
    (row) => row.getAttribute("data-component-type"),
  );
}
beforeEach(() => localStorage.clear());
afterEach(() => {
  cleanup();
  vi.restoreAllMocks();
  vi.unstubAllGlobals();
});

describe("component reordering", () => {
  it("list press dragging updates order, requests validation and persists with IDs", async () => {
    const mock = install((init) =>
      response({
        ...preview,
        component_order: JSON.parse(String(init?.body)).component_order,
      }),
    );
    render(<Workbench />);
    await waitFor(() =>
      expect(screen.getByRole("button", { name: "生成设计" })).toBeEnabled(),
    );
    listDrag();
    expect(rowOrder()).toEqual(movedOrder);
    await waitFor(() =>
      expect(screen.getByRole("button", { name: "生成设计" })).toBeEnabled(),
    );
    await waitFor(() =>
      expect(
        mock.mock.calls
          .filter(([url]) => url.includes("preview"))
          .map(([, init]) => JSON.parse(String(init?.body)).component_order),
      ).toContainEqual(movedOrder),
    );
    expect(JSON.parse(localStorage.getItem("plasmid:drag-demo")!)).toEqual({
      resistance_id: "amp",
      replication_id: "ori",
      t0_id: null,
      t1_id: null,
      component_order: movedOrder,
    });
    await act(async () =>
      fireEvent.click(screen.getByRole("button", { name: "选择 KanR" })),
    );
    expect(rowOrder()).toEqual(movedOrder);
  });

  it("ring press dragging ignores old validation and hides download", async () => {
    let resolveOld!: (response: Response) => void;
    let calls = 0;
    const mock = install(
      (init) =>
        ++calls === 1
          ? new Promise((resolve) => {
              resolveOld = resolve;
            })
          : response({
              ...preview,
              valid: false,
              component_order: JSON.parse(String(init?.body)).component_order,
            }),
      { ...context, result },
    );
    const { container } = render(<Workbench />);
    await screen.findByRole("link", { name: "下载 GenBank" });
    const svg = screen.getByRole("img", { name: "质粒环图" });
    vi.spyOn(svg, "getBoundingClientRect").mockReturnValue(
      rect(100, 50, 440, 440),
    );
    const part = container.querySelector("[data-ring-component='resistance']")!;
    expect(part).toBeInTheDocument();
    fireEvent.mouseDown(part, { button: 0, clientX: 320, clientY: 130 });
    fireEvent.mouseMove(window, { clientX: 280, clientY: 130 });
    expect(screen.getByTestId("ring-drag-ghost")).toBeInTheDocument();
    fireEvent.mouseUp(window, { clientX: 280, clientY: 130 });
    expect(rowOrder()).toEqual(movedOrder);
    expect(
      screen.queryByRole("link", { name: "下载 GenBank" }),
    ).not.toBeInTheDocument();
    await act(async () =>
      resolveOld(
        new Response(
          JSON.stringify({ ...preview, component_order: defaultOrder }),
        ),
      ),
    );
    expect(screen.getByRole("button", { name: "生成设计" })).toBeDisabled();
    expect(
      mock.mock.calls
        .filter(([url]) => url.includes("preview"))
        .map(([, init]) => JSON.parse(String(init?.body)).component_order),
    ).toContainEqual(movedOrder);
  });

  it("moves complete server blocks and site/gene coordinates immediately while preview is pending", async () => {
    let calls = 0;
    install(() =>
      ++calls === 1
        ? response({
            ...preview,
            component_order: defaultOrder,
            features: [
              ...preview.features,
              {
                label: "gene",
                type: "CDS",
                kind: "gene",
                strand: 1,
                start_bp: 67,
                end_bp: 80,
              },
            ],
          })
        : new Promise(() => {}),
    );
    const { container } = render(<Workbench />);
    await waitFor(() =>
      expect(screen.getByRole("button", { name: "生成设计" })).toBeEnabled(),
    );
    listDrag();
    expect(
      container.querySelector("[data-ring-component='resistance'] path"),
    ).toHaveAttribute("data-start-bp", "91");
    expect(
      container.querySelector("path[data-start-bp='111']"),
    ).toHaveAttribute("data-end-bp", "112");
    expect(container.querySelector("path[data-start-bp='34']")).toHaveAttribute(
      "data-end-bp",
      "38",
    );
    expect(container.querySelector("path[data-start-bp='39']")).toHaveAttribute(
      "data-end-bp",
      "44",
    );
    expect(
      container.querySelector("[data-testid='restriction-tick']"),
    ).toHaveAttribute("data-start-bp", "85");
    fireEvent.click(screen.getByRole("button", { name: "显示基因注释" }));
    expect(
      container.querySelector("path[data-start-bp='45'][data-end-bp='58']"),
    ).toBeInTheDocument();
    expect(screen.getByText(/等待验证/)).toBeInTheDocument();
    expect(screen.getByRole("button", { name: "生成设计" })).toBeDisabled();
  });

  it("list cancellation and clicks leave order unchanged and source cannot be cleared", async () => {
    const mock = install();
    const { unmount } = render(<Workbench />);
    await waitFor(() =>
      expect(screen.getByRole("button", { name: "生成设计" })).toBeEnabled(),
    );
    layoutRows();
    const row = screen.getByTestId("source-slot");
    fireEvent.mouseDown(row, { button: 0, clientX: 620, clientY: 230 });
    fireEvent.mouseMove(window, { clientX: 620, clientY: 105 });
    fireEvent.keyDown(window, { key: "Escape" });
    fireEvent.mouseUp(window, { clientX: 620, clientY: 105 });
    expect(rowOrder()).toEqual(defaultOrder);
    fireEvent.mouseDown(row, { button: 0, clientX: 620, clientY: 230 });
    fireEvent.mouseUp(window, { clientX: 620, clientY: 230 });
    expect(rowOrder()).toEqual(defaultOrder);
    expect(
      screen.queryByRole("button", { name: /清空.*表达/ }),
    ).not.toBeInTheDocument();
    fireEvent.mouseDown(row, { button: 0, clientX: 620, clientY: 230 });
    unmount();
    fireEvent.mouseUp(window, { clientX: 620, clientY: 105 });
    expect(
      mock.mock.calls.filter(([url]) => url.includes("preview")),
    ).toHaveLength(1);
  });

  it("restores saved order over older manifest selection, validates permutations, sends generation order", async () => {
    localStorage.setItem(
      "plasmid:drag-demo",
      JSON.stringify({
        resistance_id: "amp",
        replication_id: "ori",
        component_order: movedOrder,
      }),
    );
    const mock = install((init) =>
      response({
        ...preview,
        component_order: JSON.parse(String(init?.body)).component_order,
      }),
    );
    const mounted = render(<Workbench />);
    await waitFor(() =>
      expect(screen.getByRole("button", { name: "生成设计" })).toBeEnabled(),
    );
    expect(rowOrder()).toEqual(movedOrder);
    fireEvent.click(screen.getByRole("button", { name: "生成设计" }));
    await waitFor(() =>
      expect(
        mock.mock.calls
          .filter(([url]) => url.includes("generate"))
          .map(([, init]) => JSON.parse(String(init?.body)).component_order),
      ).toEqual([movedOrder]),
    );
    mounted.unmount();
    localStorage.setItem(
      "plasmid:drag-demo",
      JSON.stringify({
        resistance_id: "amp",
        replication_id: "ori",
        component_order: ["expression", "expression", "resistance"],
      }),
    );
    render(<Workbench />);
    await screen.findByRole("button", { name: "选择 AmpR" });
    expect(rowOrder()).toEqual(defaultOrder);
  });

  it("maps SVG letterboxing, rotation and zoom, and cancels on Escape and unmount", () => {
    const onReorder = vi.fn();
    const props = {
      preview,
      construct: null,
      enzymes: [],
      componentOrder: defaultOrder,
      onReorder,
    };
    const { container, unmount } = render(<PlasmidRing {...props} />);
    const svg = screen.getByRole("img", { name: "质粒环图" });
    vi.spyOn(svg, "getBoundingClientRect").mockReturnValue(
      rect(100, 50, 660, 440),
    );
    fireEvent.click(screen.getByRole("button", { name: "旋转环图" }));
    fireEvent.click(screen.getByRole("button", { name: "放大环图" }));
    const theta = -Math.PI / 2 + Math.PI * 2 * 0.95 + Math.PI / 6;
    const x = 430 + 138 * 1.1 * Math.cos(theta),
      y = 270 + 138 * 1.1 * Math.sin(theta);
    const part = container.querySelector("[data-ring-component='resistance']")!;
    expect(part).toBeInTheDocument();
    fireEvent.mouseDown(part, { button: 0, clientX: 500, clientY: 120 });
    fireEvent.mouseMove(window, { clientX: x, clientY: y });
    expect(screen.getByTestId("ring-insertion-marker")).toBeInTheDocument();
    fireEvent.mouseUp(window, { clientX: x, clientY: y });
    expect(onReorder).toHaveBeenCalledWith(movedOrder);
    expect(part).toBeInTheDocument();
    fireEvent.mouseDown(part, { button: 0, clientX: 500, clientY: 120 });
    fireEvent.mouseMove(window, { clientX: x, clientY: y });
    fireEvent.keyDown(window, { key: "Escape" });
    expect(screen.queryByTestId("ring-drag-ghost")).not.toBeInTheDocument();
    fireEvent.mouseUp(window, { clientX: x, clientY: y });
    expect(onReorder).toHaveBeenCalledTimes(1);
    expect(part).toBeInTheDocument();
    fireEvent.mouseDown(part, { button: 0, clientX: 500, clientY: 120 });
    unmount();
    fireEvent.mouseUp(window, { clientX: x, clientY: y });
    expect(onReorder).toHaveBeenCalledTimes(1);
  });

  it("keeps ring drag alive across identical context snapshots and cancels when source changes", () => {
    const onReorder = vi.fn();
    const props = {
      preview,
      construct: null,
      enzymes: [],
      componentOrder: defaultOrder,
      onReorder,
    };
    const mounted = render(<PlasmidRing {...props} />);
    const svg = screen.getByRole("img", { name: "质粒环图" });
    vi.spyOn(svg, "getBoundingClientRect").mockReturnValue(
      rect(100, 50, 440, 440),
    );
    const part = mounted.container.querySelector(
      "[data-ring-component='resistance']",
    )!;
    fireEvent.mouseDown(part, { button: 0, clientX: 320, clientY: 130 });
    fireEvent.mouseMove(window, { clientX: 280, clientY: 130 });
    mounted.rerender(
      <PlasmidRing
        {...props}
        preview={{
          ...preview,
          segments: preview.segments.map((segment) => ({ ...segment })),
        }}
        componentOrder={[...defaultOrder]}
      />,
    );
    fireEvent.mouseUp(window, { clientX: 280, clientY: 130 });
    expect(onReorder).toHaveBeenCalledWith(movedOrder);
    expect(screen.queryByTestId("ring-drag-ghost")).not.toBeInTheDocument();
    fireEvent.mouseDown(part, { button: 0, clientX: 320, clientY: 130 });
    fireEvent.mouseMove(window, { clientX: 280, clientY: 130 });
    mounted.rerender(
      <PlasmidRing
        {...props}
        preview={{ ...preview, source_fingerprint: "changed" }}
      />,
    );
    expect(screen.queryByTestId("ring-drag-ghost")).not.toBeInTheDocument();
    fireEvent.mouseUp(window, { clientX: 280, clientY: 130 });
    expect(onReorder).toHaveBeenCalledTimes(1);
  });

  it("shows selected real module lengths before selection is complete and allows drag", async () => {
    install(undefined, { ...context, selection: null });
    const { container } = render(<Workbench />);
    fireEvent.click(await screen.findByRole("button", { name: "选择 AmpR" }));
    expect(
      container.querySelector("[data-ring-component='resistance']"),
    ).toBeInTheDocument();
    expect(
      container.querySelector("[data-ring-component='expression']"),
    ).toBeInTheDocument();
    expect(screen.getByText(/等待验证/)).toBeInTheDocument();
    listDrag();
    expect(rowOrder()).toEqual(movedOrder);
    expect(screen.getByRole("button", { name: "生成设计" })).toBeDisabled();
  });
});
