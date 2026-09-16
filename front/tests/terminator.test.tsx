import { afterEach, beforeEach, expect, it, vi } from "vitest";
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
import { normalizeOrder, reorderGeometry } from "../src/componentOrder";
import type { Preview } from "../src/types";
const order = ["resistance", "replication", "t1", "expression", "t0"];
const module = (id: string, type: string, role?: string) => ({
  id,
  name: role?.toUpperCase() || id,
  type,
  role,
  sequence: "ATGC",
  length_bp: 4,
  gc_percent: 50,
  notes: [],
});
const library = {
  resistance: [module("amp", "resistance")],
  replication: [module("ori", "replication")],
  terminator: [
    module("basic_seva_t0", "terminator", "t0"),
    module("basic_seva_t1", "terminator", "t1"),
  ],
};
const context = {
  project: {
    target: "terminator-demo",
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
    length_bp: 1000,
    gc_percent: 50,
    sequence_sha256: "source",
    features: [],
  },
  restriction_enzymes: [],
  selection: null,
  result: null,
  active_job: null,
};
function install(warnings = ["未添加 T0", "未添加 T1"], valid = true) {
  const fetcher = vi.fn((url: string, init?: RequestInit) => {
    const body = JSON.parse(String(init?.body || "{}"));
    return Promise.resolve(
      new Response(
        JSON.stringify(
          url === "/api/modules"
            ? library
            : url === "/api/preview"
              ? {
                  ...body,
                  valid,
                  issues: valid
                    ? []
                    : [{ code: "enzyme", message: "酶切冲突" }],
                  warnings: [...warnings, "普通说明"],
                  terminator_warnings: warnings,
                  sequence: "ATGC",
                  length_bp: 1008,
                  gc_percent: 50,
                  sequence_sha256: "preview",
                  source_fingerprint: "source",
                  manifest_revision: 1,
                  segments: [],
                  features: [],
                  enzymes: [],
                }
              : url === "/api/generate"
                ? {
                    id: "job",
                    status: "queued",
                    stage: "queued",
                    message: "生成中",
                  }
                : context,
        ),
      ),
    );
  });
  vi.stubGlobal("fetch", fetcher);
  return fetcher;
}
async function selectRequired() {
  fireEvent.click(await screen.findByRole("button", { name: "选择 amp" }));
  fireEvent.click(screen.getByRole("button", { name: "选择 ori" }));
  await waitFor(() =>
    expect(screen.getByRole("button", { name: "生成设计" })).toBeEnabled(),
  );
}
beforeEach(() => localStorage.clear());
afterEach(() => {
  cleanup();
  vi.restoreAllMocks();
  vi.unstubAllGlobals();
});
it("keeps terminators absent until manually added and sends ordered instances", async () => {
  const fetcher = install();
  render(<Workbench />);
  await selectRequired();
  expect(screen.queryByTestId("t0-slot")).not.toBeInTheDocument();
  expect(screen.queryByTestId("t1-slot")).not.toBeInTheDocument();
  expect(
    JSON.parse(
      String(fetcher.mock.calls.find(([u]) => u === "/api/preview")?.[1]?.body),
    ),
  ).toMatchObject({
    components: [
      { component_type: "replication", module_id: "ori" },
      { component_type: "expression" },
      { component_type: "resistance", module_id: "amp" },
    ],
  });
  fireEvent.click(screen.getByRole("button", { name: "选择 T0" }));
  await waitFor(() =>
    expect(screen.getByRole("button", { name: "生成设计" })).toBeEnabled(),
  );
  expect(screen.getByTestId("t0-slot")).toHaveTextContent("T0");
  expect(screen.queryByTestId("t1-slot")).not.toBeInTheDocument();
});
it.each(["返回调整", "Escape", "backdrop"])(
  "cancels warning dialog via %s without generation or changing selections",
  async (method) => {
    const fetcher = install();
    render(<Workbench />);
    await selectRequired();
    fireEvent.click(screen.getByRole("button", { name: "生成设计" }));
    expect(screen.getByRole("alertdialog")).toHaveTextContent("未添加 T0");
    if (method === "Escape") fireEvent.keyDown(window, { key: "Escape" });
    else if (method === "backdrop")
      fireEvent.click(screen.getByTestId("terminator-backdrop"));
    else fireEvent.click(screen.getByRole("button", { name: method }));
    expect(screen.queryByRole("alertdialog")).not.toBeInTheDocument();
    expect(
      fetcher.mock.calls.filter(([u]) => u === "/api/generate"),
    ).toHaveLength(0);
    expect(screen.getByTestId("resistance-slot")).toHaveTextContent("amp");
  },
);
it("continues warnings with exact choices only once and changing selection closes dialog", async () => {
  const fetcher = install(["T1 位置不正确"]);
  render(<Workbench />);
  await selectRequired();
  fireEvent.click(screen.getByRole("button", { name: "选择 T0" }));
  await waitFor(() =>
    expect(screen.getByRole("button", { name: "生成设计" })).toBeEnabled(),
  );
  fireEvent.click(screen.getByRole("button", { name: "生成设计" }));
  fireEvent.click(screen.getByRole("button", { name: "选择 T1" }));
  expect(screen.queryByRole("alertdialog")).not.toBeInTheDocument();
  await waitFor(() =>
    expect(screen.getByRole("button", { name: "生成设计" })).toBeEnabled(),
  );
  fireEvent.click(screen.getByRole("button", { name: "生成设计" }));
  const proceed = screen.getByRole("button", { name: "继续生成" });
  fireEvent.click(proceed);
  fireEvent.click(proceed);
  await waitFor(() =>
    expect(
      fetcher.mock.calls.filter(([u]) => u === "/api/generate"),
    ).toHaveLength(1),
  );
  expect(screen.queryByRole("alertdialog")).not.toBeInTheDocument();
  expect(
    JSON.parse(
      String(
        fetcher.mock.calls.find(([u]) => u === "/api/generate")?.[1]?.body,
      ),
    ),
  ).toMatchObject({
    components: [
      { component_type: "replication", module_id: "ori" },
      { component_type: "expression" },
      { component_type: "resistance", module_id: "amp" },
      { component_type: "t0", module_id: "basic_seva_t0" },
      { component_type: "t1", module_id: "basic_seva_t1" },
    ],
  });
});
it("does not confirm generic notes and keeps hard conflicts blocking", async () => {
  const fetcher = install([]);
  const view = render(<Workbench />);
  await selectRequired();
  fireEvent.click(screen.getByRole("button", { name: "生成设计" }));
  await waitFor(() =>
    expect(fetcher.mock.calls.some(([u]) => u === "/api/generate")).toBe(true),
  );
  expect(screen.queryByRole("alertdialog")).not.toBeInTheDocument();
  view.unmount();
  install(["未添加 T0"], false);
  render(<Workbench />);
  fireEvent.click(await screen.findByRole("button", { name: "选择 amp" }));
  fireEvent.click(screen.getByRole("button", { name: "选择 ori" }));
  await screen.findByText("酶切冲突");
  expect(screen.getByRole("button", { name: "生成设计" })).toBeDisabled();
});
it("migrates legacy order around expression and respects explicit terminator segment ownership", () => {
  expect(normalizeOrder(["expression", "resistance", "replication"])).toEqual([
    "t1",
    "expression",
    "t0",
    "resistance",
    "replication",
  ]);
  const preview = {
    valid: true,
    sequence: "DNA",
    segments: [
      {
        id: "t1",
        kind: "terminator",
        component_type: "t1",
        label: "T1",
        start_bp: 1,
        end_bp: 4,
        length_bp: 4,
      },
      {
        id: "ori",
        kind: "replication",
        component_type: "replication",
        label: "ori",
        start_bp: 5,
        end_bp: 8,
        length_bp: 4,
      },
    ],
    features: [],
    length_bp: 8,
    issues: [],
    warnings: [],
    gc_percent: 50,
    sequence_sha256: "x",
    source_fingerprint: "source",
    manifest_revision: 1,
    enzymes: [],
  } as Preview;
  const next = reorderGeometry(preview, [
    "replication",
    "t1",
    "expression",
    "t0",
    "resistance",
  ]);
  expect(next.segments.find((s) => s.id === "t1")?.start_bp).toBe(5);
});

it("always labels short independent terminator sectors with external callouts and allows ring dragging", () => {
  const parts: Preview["segments"] = [
    {
      id: "ori",
      label: "ori",
      kind: "replication",
      component_type: "replication",
      start_bp: 1,
      end_bp: 990,
      length_bp: 990,
    },
    {
      id: "basic_seva_t1",
      label: "T1",
      kind: "terminator",
      component_type: "t1",
      start_bp: 991,
      end_bp: 995,
      length_bp: 5,
    },
    {
      id: "basic_seva_t0",
      label: "T0",
      kind: "terminator",
      component_type: "t0",
      start_bp: 996,
      end_bp: 1000,
      length_bp: 5,
    },
  ];
  const preview: Preview = {
    valid: true,
    issues: [],
    warnings: [],
    sequence: "",
    length_bp: 1000,
    gc_percent: 50,
    sequence_sha256: "x",
    source_fingerprint: "source",
    manifest_revision: 1,
    segments: parts,
    features: [
      {
        label: "EcoRI",
        kind: "restriction",
        type: "misc_feature",
        start_bp: 996,
        end_bp: 1000,
        strand: 1,
      },
    ],
    enzymes: [],
  };
  const reorder = vi.fn();
  const { container } = render(
    <PlasmidRing
      preview={preview}
      construct={null}
      enzymes={[]}
      onReorder={reorder}
    />,
  );
  expect(container.querySelectorAll(".terminator-callout")).toHaveLength(2);
  expect(screen.getByText("T0")).toBeVisible();
  expect(screen.getByText("T1")).toBeVisible();
  expect(screen.getByText("T0").getAttribute("y")).not.toBe(
    screen.getByText("EcoRI").getAttribute("y"),
  );
  const svg = screen.getByRole("img", { name: "质粒环图" });
  vi.spyOn(svg, "getBoundingClientRect").mockReturnValue({
    left: 0,
    top: 0,
    width: 440,
    height: 440,
  } as DOMRect);
  fireEvent.mouseDown(container.querySelector('[data-ring-component="t1"]')!, {
    button: 0,
    clientX: 215,
    clientY: 80,
  });
  fireEvent.mouseMove(window, { clientX: 320, clientY: 220 });
  fireEvent.mouseUp(window, { clientX: 320, clientY: 220 });
  expect(reorder).toHaveBeenCalledWith([
    "resistance",
    "t1",
    "replication",
    "expression",
    "t0",
  ]);
});
it("restores only explicit valid terminator references and clear all preserves source", async () => {
  localStorage.setItem(
    "plasmid:terminator-demo",
    JSON.stringify({
      resistance_id: "amp",
      replication_id: "ori",
      t0_id: "basic_seva_t0",
      t1_id: "unknown",
      component_order: ["resistance", "replication", "expression"],
    }),
  );
  install();
  render(<Workbench />);
  await waitFor(() =>
    expect(screen.getByRole("button", { name: "生成设计" })).toBeEnabled(),
  );
  expect(screen.getByTestId("t0-slot")).toHaveTextContent("T0");
  expect(screen.queryByTestId("t1-slot")).not.toBeInTheDocument();
  fireEvent.click(screen.getByRole("button", { name: "生成设计" }));
  fireEvent.click(screen.getByRole("button", { name: "清空T0" }));
  expect(screen.queryByRole("alertdialog")).not.toBeInTheDocument();
  expect(screen.getByRole("button", { name: "生成设计" })).toBeDisabled();
  fireEvent.click(screen.getByRole("button", { name: "清空" }));
  expect(
    JSON.parse(localStorage.getItem("plasmid:terminator-demo")!),
  ).toMatchObject({
    components: [
      { instance_id: "replication", component_type: "replication" },
      { instance_id: "expression", component_type: "expression" },
    ],
  });
  expect(screen.getByTestId("source-slot")).toHaveTextContent("完整表达构建");
});

it("supports manual terminator library drop and independent list reordering which closes warnings", async () => {
  const fetcher = install(["T0 位置不正确"]);
  const { container } = render(<Workbench />);
  await selectRequired();
  const values = new Map<string, string>();
  const data = {
    setData: (key: string, value: string) => values.set(key, value),
    getData: (key: string) => values.get(key) || "",
    get types() {
      return [...values.keys()];
    },
    effectAllowed: "none",
    dropEffect: "none",
  };
  const card = screen.getByRole("button", { name: "选择 T0" });
  fireEvent.dragStart(card, { dataTransfer: data });
  fireEvent.dragOver(screen.getByRole("img", { name: "质粒环图" }), {
    dataTransfer: data,
  });
  fireEvent.drop(screen.getByRole("img", { name: "质粒环图" }), {
    dataTransfer: data,
  });
  fireEvent.dragEnd(card, { dataTransfer: data });
  expect(screen.getByTestId("t0-slot")).toHaveTextContent("T0");
  expect(screen.queryByTestId("t1-slot")).not.toBeInTheDocument();
  await waitFor(() =>
    expect(screen.getByRole("button", { name: "生成设计" })).toBeEnabled(),
  );
  fireEvent.click(screen.getByRole("button", { name: "生成设计" }));
  expect(screen.getByRole("alertdialog")).toBeInTheDocument();
  ["replication", "expression", "resistance", "t0"].forEach((type, index) =>
    vi
      .spyOn(
        container.querySelector(`[data-component-type="${type}"]`)!,
        "getBoundingClientRect",
      )
      .mockReturnValue({
        left: 600,
        top: 100 + index * 60,
        width: 200,
        height: 50,
      } as DOMRect),
  );
  fireEvent.mouseDown(screen.getByTestId("t0-slot"), {
    button: 0,
    clientX: 620,
    clientY: 360,
  });
  fireEvent.mouseMove(window, { clientX: 620, clientY: 105 });
  fireEvent.mouseUp(window, { clientX: 620, clientY: 105 });
  expect(screen.queryByRole("alertdialog")).not.toBeInTheDocument();
  await waitFor(() =>
    expect(screen.getByRole("button", { name: "生成设计" })).toBeEnabled(),
  );
  expect(
    JSON.parse(
      String(
        fetcher.mock.calls
          .filter(([u]) => u === "/api/preview")
          .slice(-1)[0]?.[1]?.body,
      ),
    ),
  ).toMatchObject({
    components: [
      { component_type: "t0", module_id: "basic_seva_t0" },
      { component_type: "replication", module_id: "ori" },
      { component_type: "expression" },
      { component_type: "resistance", module_id: "amp" },
    ],
  });
});

it("ignores a late preview from before a manual terminator choice", async () => {
  const fetcher = install(["当前 T0 提示"]);
  const original = fetcher.getMockImplementation()!;
  let resolveOld!: (response: Response) => void;
  let count = 0;
  fetcher.mockImplementation((url, init) =>
    url === "/api/preview" && ++count === 1
      ? new Promise((resolve) => {
          resolveOld = resolve;
        })
      : original(url, init),
  );
  render(<Workbench />);
  fireEvent.click(await screen.findByRole("button", { name: "选择 amp" }));
  fireEvent.click(screen.getByRole("button", { name: "选择 ori" }));
  await waitFor(() => expect(count).toBe(1));
  fireEvent.click(screen.getByRole("button", { name: "选择 T0" }));
  await waitFor(() =>
    expect(screen.getByRole("button", { name: "生成设计" })).toBeEnabled(),
  );
  fireEvent.click(screen.getByRole("button", { name: "生成设计" }));
  await act(async () =>
    resolveOld(
      new Response(
        JSON.stringify({
          valid: true,
          issues: [],
          warnings: ["过期提示"],
          terminator_warnings: ["过期提示"],
          component_order: order,
          t0_id: null,
          t1_id: null,
          length_bp: 99999,
          sequence: "",
          segments: [],
          features: [],
          enzymes: [],
          gc_percent: 50,
          sequence_sha256: "old",
          source_fingerprint: "source",
          manifest_revision: 1,
        }),
      ),
    ),
  );
  expect(screen.getByRole("alertdialog")).toHaveTextContent("当前 T0 提示");
  expect(screen.getByRole("alertdialog")).not.toHaveTextContent("过期提示");
  expect(screen.queryByText("99,999 bp")).not.toBeInTheDocument();
  expect(screen.getByTestId("t0-slot")).toHaveTextContent("T0");
});
it("closes confirmation on a changed source while preserving manual choices and sends no POST", async () => {
  const fetcher = install();
  const original = fetcher.getMockImplementation()!;
  let changed = false;
  fetcher.mockImplementation((url, init) =>
    url === "/api/context" && changed
      ? Promise.resolve(
          new Response(
            JSON.stringify({
              ...context,
              project: {
                ...context.project,
                source_fingerprint: "new-source",
                manifest_revision: 2,
              },
            }),
          ),
        )
      : original(url, init),
  );
  render(<Workbench />);
  await selectRequired();
  fireEvent.click(screen.getByRole("button", { name: "选择 T0" }));
  await waitFor(() =>
    expect(screen.getByRole("button", { name: "生成设计" })).toBeEnabled(),
  );
  fireEvent.click(screen.getByRole("button", { name: "生成设计" }));
  expect(screen.getByRole("alertdialog")).toBeInTheDocument();
  changed = true;
  await act(async () =>
    fireEvent.click(screen.getByRole("button", { name: "刷新源构建" })),
  );
  expect(screen.queryByRole("alertdialog")).not.toBeInTheDocument();
  expect(screen.getByTestId("t0-slot")).toHaveTextContent("T0");
  expect(
    fetcher.mock.calls.filter(([u]) => u === "/api/generate"),
  ).toHaveLength(0);
  expect(
    JSON.parse(
      String(
        fetcher.mock.calls
          .filter(([u]) => u === "/api/preview")
          .slice(-1)[0]?.[1]?.body,
      ),
    ),
  ).toMatchObject({
    source_fingerprint: "new-source",
    expected_revision: 2,
    components: expect.arrayContaining([
      {
        instance_id: expect.any(String),
        component_type: "t0",
        module_id: "basic_seva_t0",
      },
    ]),
  });
});
