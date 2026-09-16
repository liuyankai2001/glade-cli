import { afterEach, beforeEach, expect, it, vi } from "vitest";
import {
  cleanup,
  fireEvent,
  render,
  screen,
  waitFor,
  within,
} from "@testing-library/react";
import "@testing-library/jest-dom/vitest";
import { Workbench } from "../src/Workbench";

const part = (id: string, type: string, role?: string) => ({
  id,
  name: id,
  type,
  role,
  sequence: "ATGC",
  length_bp: 4,
  gc_percent: 50,
  notes: [],
});
const library = {
  resistance: [part("amp", "resistance"), part("kan", "resistance")],
  replication: [part("ori", "replication")],
  terminator: [part("T0", "terminator", "t0"), part("T1", "terminator", "t1")],
};
const context = {
  project: {
    target: "copies",
    name: "Copies",
    host: "E. coli",
    manifest_revision: 1,
    source_fingerprint: "source",
  },
  ready: true,
  issues: [],
  selection: null,
  result: null,
  active_job: null,
  construct: {
    id: 1,
    name: "Source",
    length_bp: 100,
    gc_percent: 50,
    sequence_sha256: "source",
    features: [],
  },
  restriction_enzymes: [],
};
function install() {
  vi.stubGlobal(
    "fetch",
    vi.fn((url: string, init?: RequestInit) => {
      const body = JSON.parse(String(init?.body || "{}"));
      let cursor = 1;
      const preview = {
        ...body,
        valid: true,
        issues: [],
        warnings: [],
        terminator_warnings: [],
        sequence: "",
        sequence_sha256: "preview",
        source_fingerprint: "source",
        manifest_revision: 1,
        gc_percent: 50,
        features: [],
        enzymes: [],
        segments: (body.components || []).map(
          (c: {
            instance_id: string;
            component_type: string;
            module_id?: string;
          }) => {
            const size = c.component_type === "expression" ? 100 : 4;
            const start = cursor;
            cursor += size;
            return {
              id: c.module_id || "expression",
              instance_id: c.instance_id,
              component_type: c.component_type,
              label: c.module_id || "Source",
              kind: ["t0", "t1"].includes(c.component_type)
                ? "terminator"
                : c.component_type,
              start_bp: start,
              end_bp: cursor - 1,
              length_bp: size,
            };
          },
        ),
        length_bp: 0,
      };
      preview.length_bp = cursor - 1;
      return Promise.resolve(
        new Response(
          JSON.stringify(
            url === "/api/modules"
              ? library
              : url === "/api/preview"
                ? preview
                : context,
          ),
        ),
      );
    }),
  );
}
beforeEach(() => {
  localStorage.clear();
  install();
});
afterEach(() => {
  cleanup();
  vi.restoreAllMocks();
  vi.unstubAllGlobals();
});

it("adds every click as a separate copy and removes only the chosen copy", async () => {
  const { container } = render(<Workbench />);
  const amp = await screen.findByRole("button", { name: "选择 amp" });
  fireEvent.click(amp);
  fireEvent.click(amp);
  fireEvent.click(amp);
  const t0 = screen.getByRole("button", { name: "选择 T0" });
  fireEvent.click(t0);
  fireEvent.click(t0);
  fireEvent.click(t0);
  expect(screen.getAllByTestId("resistance-slot")).toHaveLength(3);
  expect(screen.getAllByTestId("t0-slot")).toHaveLength(3);
  expect(
    container.querySelectorAll('[data-ring-component="resistance"]'),
  ).toHaveLength(3);
  expect(container.querySelectorAll('[data-ring-component="t0"]')).toHaveLength(
    3,
  );
  const rows = screen.getAllByTestId("resistance-slot");
  const ids = rows.map((row) =>
    row.closest(".component-row")!.getAttribute("data-component-instance"),
  );
  expect(new Set(ids).size).toBe(3);
  fireEvent.click(within(rows[1]).getByRole("button"));
  expect(screen.getAllByTestId("resistance-slot")).toHaveLength(2);
  expect(
    screen
      .getAllByTestId("resistance-slot")
      .map((row) =>
        row.closest(".component-row")!.getAttribute("data-component-instance"),
      ),
  ).toEqual([ids[0], ids[2]]);
  expect(screen.getAllByTestId("t0-slot")).toHaveLength(3);
  expect(
    container.querySelectorAll('[data-ring-component="resistance"]'),
  ).toHaveLength(2);
});

it("keeps different markers and all copies after a reload, with one replication module", async () => {
  const first = render(<Workbench />);
  fireEvent.click(await screen.findByRole("button", { name: "选择 amp" }));
  fireEvent.click(screen.getByRole("button", { name: "选择 amp" }));
  fireEvent.click(screen.getByRole("button", { name: "选择 kan" }));
  fireEvent.click(screen.getByRole("button", { name: "选择 ori" }));
  fireEvent.click(screen.getByRole("button", { name: "选择 ori" }));
  await waitFor(() =>
    expect(screen.getByRole("button", { name: "生成设计" })).toBeEnabled(),
  );
  const saved = JSON.parse(localStorage.getItem("plasmid:copies")!);
  expect(
    saved.components
      .filter(
        (c: { component_type: string }) => c.component_type === "resistance",
      )
      .map((c: { module_id: string }) => c.module_id),
  ).toEqual(["amp", "amp", "kan"]);
  first.unmount();
  render(<Workbench />);
  await waitFor(() =>
    expect(screen.getAllByTestId("resistance-slot")).toHaveLength(3),
  );
  expect(screen.getAllByTestId("replication-slot")).toHaveLength(1);
  expect(
    JSON.parse(localStorage.getItem("plasmid:copies")!).components,
  ).toEqual(saved.components);
});

it("moves just one duplicate in the list and keeps the ring in the same order", async () => {
  const { container } = render(<Workbench />);
  const amp = await screen.findByRole("button", { name: "选择 amp" });
  fireEvent.click(amp);
  fireEvent.click(amp);
  const rows = Array.from(
    container.querySelectorAll<HTMLElement>(".component-row"),
  );
  const ids = rows.map((row) => row.dataset.componentInstance);
  rows.forEach((row, index) =>
    vi.spyOn(row, "getBoundingClientRect").mockReturnValue({
      x: 0,
      y: index * 50,
      left: 0,
      top: index * 50,
      width: 200,
      height: 40,
      right: 200,
      bottom: index * 50 + 40,
      toJSON() {},
    }),
  );
  fireEvent.mouseDown(rows[rows.length - 1], {
    button: 0,
    clientX: 20,
    clientY: 160,
  });
  fireEvent.mouseMove(window, { clientX: 20, clientY: -20 });
  fireEvent.mouseUp(window, { clientX: 20, clientY: -20 });
  expect(
    Array.from(container.querySelectorAll<HTMLElement>(".component-row")).map(
      (row) => row.dataset.componentInstance,
    ),
  ).toEqual([ids[3], ids[0], ids[1], ids[2]]);
  expect(
    Array.from(
      container.querySelectorAll<HTMLElement>("[data-ring-instance]"),
    ).map((row) => row.dataset.ringInstance),
  ).toEqual([ids[3], "expression", ids[2]]);
});

it("colors the ring drag feedback by module kind and moves only the chosen copy", async () => {
  const { container } = render(<Workbench />);
  const amp = await screen.findByRole("button", { name: "选择 amp" });
  fireEvent.click(amp);
  fireEvent.click(amp);
  const parts = container.querySelectorAll(
    '[data-ring-component="resistance"]',
  );
  const ids = Array.from(parts).map((part) =>
    part.getAttribute("data-ring-instance"),
  );
  const svg = screen.getByRole("img", { name: "质粒环图" });
  vi.spyOn(svg, "getBoundingClientRect").mockReturnValue({
    x: 0,
    y: 0,
    left: 0,
    top: 0,
    width: 440,
    height: 440,
    right: 440,
    bottom: 440,
    toJSON() {},
  });
  fireEvent.mouseDown(parts[1], { button: 0, clientX: 210, clientY: 80 });
  fireEvent.mouseMove(window, { clientX: 320, clientY: 220 });
  expect(screen.getByTestId("ring-drag-ghost")).toHaveAttribute(
    "fill",
    "#3fb950",
  );
  fireEvent.mouseUp(window, { clientX: 320, clientY: 220 });
  expect(
    Array.from(container.querySelectorAll<HTMLElement>(".component-row")).map(
      (row) => row.dataset.componentInstance,
    ),
  ).toEqual(["replication", ids[1], "expression", ids[0]]);
  expect(screen.getAllByTestId("resistance-slot")).toHaveLength(2);
});
