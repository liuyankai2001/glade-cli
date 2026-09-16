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
import type { Preview } from "../src/types";
import { PlasmidRing, Workbench } from "../src/Workbench";

const context = {
  project: {
    target: "demo",
    name: "Demo construct",
    host: "E. coli",
    manifest_revision: 3,
    source_fingerprint: "source-a",
  },
  ready: true,
  issues: [],
  construct: {
    id: 1,
    name: "Source",
    length_bp: 2000,
    gc_percent: 50,
    sequence_sha256: "abc",
    features: [],
  },
  restriction_enzymes: [
    { name: "BamHI", site: "GGATCC", role: "left" as const, overhang: "GATC" },
    { name: "EcoRI", site: "GAATTC", role: "right" as const, overhang: "AATT" },
  ],
  selection: null,
  result: null,
  active_job: null,
};
const modules = {
  resistance: [
    {
      id: "amp",
      name: "AmpR",
      type: "resistance",
      sequence: "ATGC",
      length_bp: 861,
      gc_percent: 53,
      antibiotic: "Ampicillin",
      resistance_gene: "bla",
      notes: ["37°C"],
    },
  ],
  replication: [
    {
      id: "puc",
      name: "pUC ori",
      type: "replication",
      sequence: "ATGC",
      length_bp: 520,
      gc_percent: 51,
      host_range: "E. coli",
      copy_number: "high",
      notes: ["高拷贝"],
    },
  ],
};
const preview: Preview = {
  valid: true,
  issues: [],
  warnings: [],
  sequence: "ATGC",
  length_bp: 3381,
  gc_percent: 51,
  sequence_sha256: "next",
  segments: [
    {
      id: "amp",
      label: "AmpR",
      kind: "resistance",
      start_bp: 1,
      end_bp: 861,
      length_bp: 861,
    },
    {
      id: "puc",
      label: "pUC ori",
      kind: "replication",
      start_bp: 862,
      end_bp: 1381,
      length_bp: 520,
    },
  ],
  features: [],
  enzymes: context.restriction_enzymes,
  source_fingerprint: "source-a",
  manifest_revision: 3,
};

function response(body: unknown, ok = true, status = 200) {
  return Promise.resolve(
    new Response(JSON.stringify(body), {
      status,
      headers: { "Content-Type": "application/json" },
    }),
  );
}
function transfer() {
  const data = new Map<string, string>();
  return {
    setData: (type: string, value: string) => data.set(type, value),
    getData: (type: string) => data.get(type) || "",
    get types() {
      return [...data.keys()];
    },
    effectAllowed: "none",
    dropEffect: "none",
  };
}
beforeEach(() => {
  localStorage.clear();
});
afterEach(() => {
  cleanup();
  vi.restoreAllMocks();
});

describe("Workbench", () => {
  it.each(["ring", "list"])(
    "accepts library dragging onto the whole %s",
    async (target) => {
      const mock = vi.fn((url: string, _init?: RequestInit) =>
        response(
          url.includes("modules")
            ? modules
            : url.includes("preview")
              ? preview
              : context,
        ),
      );
      vi.stubGlobal("fetch", mock);
      render(<Workbench />);
      const amp = await screen.findByRole("button", { name: "选择 AmpR" });
      const destination =
        target === "ring"
          ? screen.getByRole("img", { name: "质粒环图" })
          : screen.getByText("📋 组装列表");
      const ampData = transfer();
      fireEvent.dragStart(amp, { dataTransfer: ampData });
      expect(ampData.effectAllowed).toBe("copy");
      fireEvent.dragOver(destination, { dataTransfer: ampData });
      expect(ampData.dropEffect).toBe("copy");
      fireEvent.drop(destination, { dataTransfer: ampData });
      fireEvent.dragEnd(amp, { dataTransfer: ampData });
      expect(screen.getByTestId("resistance-slot")).toHaveTextContent("AmpR");
      expect(screen.getByTestId("source-slot")).toHaveTextContent(
        "完整表达构建",
      );
      const repData = transfer();
      const rep = screen.getByRole("button", { name: "选择 pUC ori" });
      fireEvent.dragStart(rep, { dataTransfer: repData });
      fireEvent.dragOver(destination, { dataTransfer: repData });
      fireEvent.drop(destination, { dataTransfer: repData });
      fireEvent.dragEnd(rep, { dataTransfer: repData });
      expect(screen.getByTestId("replication-slot")).toHaveTextContent(
        "pUC ori",
      );
      await waitFor(() =>
        expect(screen.getByRole("button", { name: "生成设计" })).toBeEnabled(),
      );
      const request = mock.mock.calls.find(([url]) => url.includes("preview"));
      expect(JSON.parse(String(request?.[1]?.body))).toEqual({
        component_order: [
          "resistance",
          "replication",
          "t1",
          "expression",
          "t0",
        ],
        t0_id: null,
        t1_id: null,
        resistance_id: "amp",
        replication_id: "puc",
        expected_revision: 3,
        source_fingerprint: "source-a",
      });
    },
  );
  it("replaces a selected module by dropping another of the same type onto the ring", async () => {
    vi.stubGlobal(
      "fetch",
      vi.fn((url: string) =>
        response(
          url.includes("modules")
            ? {
                ...modules,
                resistance: [
                  ...modules.resistance,
                  { ...modules.resistance[0], id: "kan", name: "KanR" },
                ],
              }
            : url.includes("preview")
              ? preview
              : context,
        ),
      ),
    );
    render(<Workbench />);
    fireEvent.click(await screen.findByRole("button", { name: "选择 AmpR" }));
    await act(async () =>
      fireEvent.click(screen.getByRole("button", { name: "选择 pUC ori" })),
    );
    const data = transfer();
    const kan = screen.getByRole("button", { name: "选择 KanR" });
    fireEvent.dragStart(kan, { dataTransfer: data });
    await act(async () =>
      fireEvent.drop(screen.getByRole("img", { name: "质粒环图" }), {
        dataTransfer: data,
      }),
    );
    expect(screen.getByTestId("resistance-slot")).toHaveTextContent("KanR");
    expect(screen.getByTestId("resistance-slot")).not.toHaveTextContent("AmpR");
    expect(screen.getByTestId("replication-slot")).toHaveTextContent("pUC ori");
  });
  it("ignores foreign drag data and rejects dropping into a slot of the wrong type", async () => {
    vi.stubGlobal(
      "fetch",
      vi.fn((url: string) =>
        response(url.includes("modules") ? modules : context),
      ),
    );
    render(<Workbench />);
    const amp = await screen.findByRole("button", { name: "选择 AmpR" });
    const data = transfer();
    fireEvent.dragStart(amp, { dataTransfer: data });
    fireEvent.drop(screen.getByTestId("replication-slot"), {
      dataTransfer: data,
    });
    const foreign = transfer();
    foreign.setData("text/plain", "unknown-module");
    foreign.setData(
      "application/json",
      JSON.stringify({ id: "amp", sequence: "AAAA" }),
    );
    fireEvent.drop(screen.getByRole("img", { name: "质粒环图" }), {
      dataTransfer: foreign,
    });
    expect(screen.getByTestId("resistance-slot")).not.toHaveTextContent("AmpR");
    expect(screen.getByTestId("replication-slot")).not.toHaveTextContent(
      "pUC ori",
    );
    expect(screen.getByRole("button", { name: "生成设计" })).toBeDisabled();
  });
  it("replaces the fixed resistance slot when another resistance is chosen", async () => {
    const second = {
      ...modules,
      resistance: [
        ...modules.resistance,
        { ...modules.resistance[0], id: "kan", name: "KanR" },
      ],
    };
    vi.stubGlobal(
      "fetch",
      vi
        .fn()
        .mockImplementation((url: string) =>
          response(url.includes("modules") ? second : context),
        ),
    );
    render(<Workbench />);
    await screen.findByText("AmpR");
    fireEvent.click(screen.getByRole("button", { name: /AmpR/ }));
    await act(async () =>
      fireEvent.click(screen.getByRole("button", { name: /KanR/ })),
    );
    expect(screen.getByTestId("resistance-slot")).toHaveTextContent("KanR");
    expect(screen.getByTestId("resistance-slot")).not.toHaveTextContent("AmpR");
  });

  it("does not generate before a valid preview exists", async () => {
    vi.stubGlobal(
      "fetch",
      vi
        .fn()
        .mockImplementation((url: string) =>
          response(url.includes("modules") ? modules : context),
        ),
    );
    render(<Workbench />);
    await screen.findByText("AmpR");
    expect(screen.getByRole("button", { name: /生成设计/ })).toBeDisabled();
  });

  it("ignores a stale preview response after selections change", async () => {
    let resolveOld!: (value: Response) => void;
    const oldPreview = new Promise<Response>((resolve) => {
      resolveOld = resolve;
    });
    const second = {
      ...modules,
      resistance: [
        ...modules.resistance,
        { ...modules.resistance[0], id: "kan", name: "KanR" },
      ],
    };
    const fetchMock = vi.fn().mockImplementation((url: string) => {
      if (url.includes("modules")) return response(second);
      if (url.includes("preview"))
        return fetchMock.mock.calls.filter(([u]) =>
          String(u).includes("preview"),
        ).length === 1
          ? oldPreview
          : response({ ...preview, length_bp: 4000 });
      return response(context);
    });
    vi.stubGlobal("fetch", fetchMock);
    render(<Workbench />);
    await screen.findByText("AmpR");
    fireEvent.click(screen.getByRole("button", { name: /AmpR/ }));
    fireEvent.click(screen.getByRole("button", { name: /pUC ori/ }));
    await act(async () =>
      fireEvent.click(screen.getByRole("button", { name: /KanR/ })),
    );
    await waitFor(() =>
      expect(screen.getByText(/4,000 bp/)).toBeInTheDocument(),
    );
    resolveOld(await response({ ...preview, length_bp: 3381 }));
    await new Promise((resolve) => setTimeout(resolve, 0));
    expect(screen.queryByText(/3,381 bp/)).not.toBeInTheDocument();
  });

  it("shows API failures and exposes actual downloaded files after a successful job", async () => {
    const fetchMock = vi
      .fn()
      .mockImplementation((url: string, init?: RequestInit) => {
        if (url.includes("modules")) return response(modules);
        if (url.includes("preview")) return response(preview);
        if (url.includes("generate"))
          return response(
            {
              id: "job-1",
              status: "succeeded",
              stage: "done",
              message: "完成",
              result: {
                id: "g1",
                resistance_id: "amp",
                replication_id: "puc",
                sequence_sha256: "next",
                length_bp: 3381,
                gc_percent: 51,
                manifest_revision: 3,
                source_fingerprint: "source-a",
                files: [
                  {
                    id: "gb",
                    label: "GenBank",
                    filename: "design.gb",
                    format: "gb",
                    url: "/api/files/g1/gb",
                  },
                ],
                warnings: [],
              },
            },
            true,
            202,
          );
        return response({
          ...context,
          project: { ...context.project, manifest_revision: 2 },
        });
      });
    vi.stubGlobal("fetch", fetchMock);
    render(<Workbench />);
    await screen.findByText("AmpR");
    fireEvent.click(screen.getByRole("button", { name: /AmpR/ }));
    fireEvent.click(screen.getByRole("button", { name: /pUC ori/ }));
    await screen.findByText(/3,381 bp/);
    fireEvent.click(screen.getByRole("button", { name: /生成设计/ }));
    expect(
      await screen.findByRole("link", { name: /下载 GenBank/ }),
    ).toHaveAttribute("href", "/api/files/g1/gb");
  });

  it("renders a failed preview API message", async () => {
    vi.stubGlobal(
      "fetch",
      vi.fn().mockImplementation((url: string) => {
        if (url.includes("modules")) return response(modules);
        if (url.includes("preview"))
          return response(
            { error: { code: "invalid_modules", message: "模块不兼容" } },
            false,
            422,
          );
        return response(context);
      }),
    );
    render(<Workbench />);
    await screen.findByText("AmpR");
    fireEvent.click(screen.getByRole("button", { name: /AmpR/ }));
    fireEvent.click(screen.getByRole("button", { name: /pUC ori/ }));
    expect(await screen.findByText(/预览失败：模块不兼容/)).toBeInTheDocument();
  });

  it("hides a completed result after a dirty selection change", async () => {
    const second = {
      ...modules,
      resistance: [
        ...modules.resistance,
        { ...modules.resistance[0], id: "kan", name: "KanR" },
      ],
    };
    const completed = {
      ...context,
      result: {
        id: "g1",
        resistance_id: "amp",
        replication_id: "puc",
        sequence_sha256: "next",
        length_bp: 3381,
        gc_percent: 51,
        manifest_revision: 3,
        source_fingerprint: "source-a",
        files: [
          {
            id: "gb",
            label: "GenBank",
            filename: "design.gb",
            format: "gb",
            url: "/api/files/g1/gb",
          },
        ],
        warnings: [],
      },
    };
    vi.stubGlobal(
      "fetch",
      vi
        .fn()
        .mockImplementation((url: string) =>
          response(
            url.includes("modules")
              ? second
              : url.includes("preview")
                ? preview
                : completed,
          ),
        ),
    );
    render(<Workbench />);
    await screen.findByText("AmpR");
    fireEvent.click(screen.getByRole("button", { name: /AmpR/ }));
    fireEvent.click(screen.getByRole("button", { name: /pUC ori/ }));
    expect(
      await screen.findByRole("link", { name: /下载 GenBank/ }),
    ).toBeInTheDocument();
    await act(async () =>
      fireEvent.click(screen.getByRole("button", { name: /KanR/ })),
    );
    expect(
      screen.queryByRole("link", { name: /下载 GenBank/ }),
    ).not.toBeInTheDocument();
  });

  it("draws full and short macro arcs from their supplied coordinates and renders source features", () => {
    const { container, getByText } = render(
      <PlasmidRing
        preview={{
          ...preview,
          length_bp: 1000,
          segments: [
            {
              id: "full",
              label: "完整环",
              kind: "gene",
              start_bp: 1,
              end_bp: 1000,
              length_bp: 1000,
            },
            {
              id: "short",
              label: "短片段",
              kind: "linker",
              start_bp: 500,
              end_bp: 510,
              length_bp: 11,
            },
          ],
          features: [
            {
              label: "source CDS",
              type: "CDS",
              kind: "gene",
              start_bp: 100,
              end_bp: 300,
              strand: 1,
            },
          ],
          enzymes: [],
        }}
        construct={null}
        enzymes={[]}
      />,
    );
    expect(getByText("完整环")).toBeInTheDocument();
    fireEvent.click(screen.getByRole("button", { name: "显示基因注释" }));
    expect(
      container.querySelector("path[data-start-bp='100']"),
    ).toHaveAttribute("data-end-bp", "300");
    expect(
      container.querySelector("[data-testid=full-segment]"),
    ).toHaveAttribute("r", "138");
    expect(container.querySelectorAll("path")).toHaveLength(2);
    expect(
      container.querySelector("[data-testid=strand-arrow]"),
    ).toHaveAttribute("data-strand", "1");
  });

  it("keeps a generated download when its source and selected IDs match despite a newer result revision", async () => {
    const result = {
      id: "g2",
      resistance_id: "amp",
      replication_id: "puc",
      sequence_sha256: "next",
      length_bp: 3381,
      gc_percent: 51,
      manifest_revision: 2,
      source_fingerprint: "source-a",
      files: [
        {
          id: "gb",
          label: "GenBank",
          filename: "design.gb",
          format: "gb",
          url: "/api/files/g2/gb",
        },
      ],
      warnings: [],
    };
    vi.stubGlobal(
      "fetch",
      vi.fn().mockImplementation((url: string) =>
        response(
          url.includes("modules")
            ? modules
            : url.includes("preview")
              ? preview
              : url.includes("generate")
                ? {
                    id: "j2",
                    status: "succeeded",
                    stage: "done",
                    message: "完成",
                    result,
                  }
                : {
                    ...context,
                    project: { ...context.project, manifest_revision: 1 },
                  },
        ),
      ),
    );
    render(<Workbench />);
    await screen.findByText("AmpR");
    fireEvent.click(screen.getByRole("button", { name: /AmpR/ }));
    fireEvent.click(screen.getByRole("button", { name: /pUC ori/ }));
    await screen.findByText(/3,381 bp/);
    fireEvent.click(screen.getByRole("button", { name: /生成设计/ }));
    expect(
      await screen.findByRole("link", { name: /下载 GenBank/ }),
    ).toHaveAttribute("href", "/api/files/g2/gb");
  });

  it("does not request preview until a source construct is ready", async () => {
    vi.stubGlobal(
      "fetch",
      vi
        .fn()
        .mockImplementation((url: string) =>
          response(
            url.includes("modules") ? modules : { ...context, ready: false },
          ),
        ),
    );
    render(<Workbench />);
    await screen.findByText("AmpR");
    fireEvent.click(screen.getByRole("button", { name: /AmpR/ }));
    fireEvent.click(screen.getByRole("button", { name: /pUC ori/ }));
    expect(
      (fetch as ReturnType<typeof vi.fn>).mock.calls.some(([url]) =>
        String(url).includes("/preview"),
      ),
    ).toBe(false);
  });
  it("recovers a generate conflict by refreshing context and preview", async () => {
    let revision = 3;
    const mock = vi.fn((url: string) => {
      if (url.includes("modules")) return response(modules);
      if (url.includes("preview"))
        return response({ ...preview, manifest_revision: revision });
      if (url.includes("generate")) {
        revision = 4;
        return response({ error: { message: "源已更新" } }, false, 409);
      }
      return response({
        ...context,
        project: { ...context.project, manifest_revision: revision },
      });
    });
    vi.stubGlobal("fetch", mock);
    render(<Workbench />);
    await screen.findByText("AmpR");
    fireEvent.click(screen.getByRole("button", { name: "选择 AmpR" }));
    fireEvent.click(screen.getByRole("button", { name: "选择 pUC ori" }));
    await waitFor(() =>
      expect(screen.getByRole("button", { name: "生成设计" })).toBeEnabled(),
    );
    fireEvent.click(screen.getByRole("button", { name: "生成设计" }));
    await waitFor(() =>
      expect(
        mock.mock.calls.filter(([url]) => url.includes("preview")).length,
      ).toBeGreaterThan(1),
    );
    await waitFor(() =>
      expect(screen.getByRole("button", { name: "生成设计" })).toBeEnabled(),
    );
  });
  it("guards pending submission, exposes read-only DNA, and supports clearing", async () => {
    const mock = vi.fn((url: string) =>
      url.includes("generate")
        ? new Promise<Response>(() => {})
        : response(
            url.includes("modules")
              ? modules
              : url.includes("preview")
                ? preview
                : context,
          ),
    );
    vi.stubGlobal("fetch", mock);
    render(<Workbench />);
    await screen.findByText("AmpR");
    fireEvent.click(screen.getByRole("button", { name: "选择 AmpR" }));
    fireEvent.click(screen.getByRole("button", { name: "选择 pUC ori" }));
    await waitFor(() =>
      expect(screen.getByRole("button", { name: "生成设计" })).toBeEnabled(),
    );
    expect(screen.getByLabelText("AmpR 只读 DNA")).toHaveAttribute("readonly");
    fireEvent.click(screen.getByRole("button", { name: "生成设计" }));
    fireEvent.click(screen.getByRole("button", { name: "生成设计" }));
    expect(
      mock.mock.calls.filter(([url]) => url.includes("generate")),
    ).toHaveLength(1);
    expect(screen.getByRole("button", { name: "生成设计" })).toBeDisabled();
    fireEvent.click(screen.getByRole("button", { name: "清空抗性标记" }));
    expect(screen.getByTestId("resistance-slot")).not.toHaveTextContent("AmpR");
    expect(screen.getByTestId("source-slot")).toHaveTextContent("完整表达构建");
  });
  it("discards cached IDs that do not exist in the current library", async () => {
    localStorage.setItem(
      "plasmid:demo",
      JSON.stringify({ resistance_id: "missing", replication_id: "puc" }),
    );
    const mock = vi.fn((url: string) =>
      response(url.includes("modules") ? modules : context),
    );
    vi.stubGlobal("fetch", mock);
    render(<Workbench />);
    await screen.findByText("AmpR");
    expect(screen.getByTestId("resistance-slot")).not.toHaveTextContent("AmpR");
    expect(screen.queryByLabelText("AmpR 只读 DNA")).not.toBeInTheDocument();
    expect(mock.mock.calls.some(([url]) => url.includes("preview"))).toBe(
      false,
    );
  });
  it("uses actual restriction coordinates and omits duplicate module inner arcs", () => {
    const { container } = render(
      <PlasmidRing
        preview={{
          ...preview,
          features: [
            {
              label: "AmpR",
              type: "CDS",
              kind: "resistance",
              start_bp: 1,
              end_bp: 861,
              strand: 1,
            },
            {
              label: "ori",
              type: "rep_origin",
              kind: "replication",
              start_bp: 862,
              end_bp: 1381,
              strand: 0,
            },
            {
              label: "EcoRI",
              type: "right_retained_site",
              kind: "restriction",
              start_bp: 2000,
              end_bp: 2005,
              strand: 1,
            },
          ],
        }}
        construct={null}
        enzymes={context.restriction_enzymes as any}
      />,
    );
    expect(
      container.querySelector("[data-testid=restriction-tick]"),
    ).toHaveAttribute("data-start-bp", "2000");
    expect(
      container.querySelectorAll("[data-testid=restriction-tick]"),
    ).toHaveLength(1);
  });
  it("refreshes the preview revision after committing a job while preserving chosen IDs", async () => {
    let revision = 1;
    const committed = {
      id: "g",
      resistance_id: "amp",
      replication_id: "puc",
      source_fingerprint: "source-a",
      manifest_revision: 2,
      files: [
        {
          id: "gb",
          label: "GenBank",
          url: "/api/files/g/gb",
          filename: "design.gb",
        },
      ],
    };
    const mock = vi.fn((url: string, init?: RequestInit) => {
      if (url.includes("modules")) return response(modules);
      if (url.includes("preview"))
        return response({ ...preview, manifest_revision: revision });
      if (url.includes("generate")) {
        revision = 2;
        return response({
          id: "j",
          status: "succeeded",
          stage: "done",
          message: "完成",
          result: committed,
        });
      }
      return response({
        ...context,
        result: revision === 2 ? committed : null,
        project: { ...context.project, manifest_revision: revision },
      });
    });
    vi.stubGlobal("fetch", mock);
    render(<Workbench />);
    await screen.findByText("AmpR");
    fireEvent.click(screen.getByRole("button", { name: "选择 AmpR" }));
    fireEvent.click(screen.getByRole("button", { name: "选择 pUC ori" }));
    await waitFor(() =>
      expect(screen.getByRole("button", { name: "生成设计" })).toBeEnabled(),
    );
    fireEvent.click(screen.getByRole("button", { name: "生成设计" }));
    await screen.findByRole("link", { name: "下载 GenBank" });
    await waitFor(() =>
      expect(
        mock.mock.calls
          .filter(([url]) => url.includes("preview"))
          .some(
            ([, init]) =>
              JSON.parse(String(init?.body)).expected_revision === 2,
          ),
      ).toBe(true),
    );
    expect(screen.getByTestId("resistance-slot")).toHaveTextContent("AmpR");
  });
  it("polls changed source context and disables generation when it becomes unready", async () => {
    let reads = 0;
    vi.stubGlobal(
      "fetch",
      vi.fn((url: string) =>
        response(
          url.includes("modules")
            ? modules
            : url.includes("preview")
              ? preview
              : { ...context, ready: ++reads === 1 },
        ),
      ),
    );
    render(<Workbench />);
    await screen.findByText("AmpR");
    fireEvent.click(screen.getByRole("button", { name: "选择 AmpR" }));
    fireEvent.click(screen.getByRole("button", { name: "选择 pUC ori" }));
    await waitFor(() =>
      expect(screen.getByRole("button", { name: "生成设计" })).toBeEnabled(),
    );
    await waitFor(
      () =>
        expect(screen.getByRole("button", { name: "生成设计" })).toBeDisabled(),
      { timeout: 6000 },
    );
  }, 7000);
  it.each(["replacement", "invalid"])(
    "uses authoritative context after a completed job: %s",
    async (mode) => {
      const resultA = {
        id: "a",
        resistance_id: "amp",
        replication_id: "puc",
        source_fingerprint: "source-a",
        manifest_revision: 2,
        files: [
          {
            id: "gb",
            label: "GenBank",
            url: "/api/files/a/gb",
            filename: "a.gb",
          },
        ],
      };
      let authoritative: any = null;
      let revision = 1;
      let invalid = false;
      vi.stubGlobal(
        "fetch",
        vi.fn((url: string) => {
          if (url.includes("modules")) return response(modules);
          if (url.includes("preview"))
            return response({ ...preview, manifest_revision: revision });
          if (url.includes("generate")) {
            revision = 2;
            authoritative = resultA;
            return response({
              id: "job-a",
              status: "succeeded",
              stage: "done",
              message: "完成",
              result: resultA,
            });
          }
          return response({
            ...context,
            project: { ...context.project, manifest_revision: revision },
            result: authoritative,
            issues: invalid
              ? [{ code: "artifact_invalid", message: "产物校验失败" }]
              : [],
          });
        }),
      );
      render(<Workbench />);
      await screen.findByText("AmpR");
      await act(async () => {
        fireEvent.click(screen.getByRole("button", { name: "选择 AmpR" }));
        fireEvent.click(screen.getByRole("button", { name: "选择 pUC ori" }));
      });
      await act(async () =>
        fireEvent.click(screen.getByRole("button", { name: "生成设计" })),
      );
      expect(
        screen.getByRole("link", { name: "下载 GenBank" }),
      ).toHaveAttribute("href", "/api/files/a/gb");
      authoritative =
        mode === "replacement"
          ? {
              ...resultA,
              id: "b",
              files: [
                {
                  id: "gb",
                  label: "GenBank",
                  url: "/api/files/b/gb",
                  filename: "b.gb",
                },
              ],
            }
          : null;
      invalid = mode === "invalid";
      await act(async () =>
        fireEvent.click(screen.getByRole("button", { name: "刷新源构建" })),
      );
      if (mode === "replacement")
        expect(
          screen.getByRole("link", { name: "下载 GenBank" }),
        ).toHaveAttribute("href", "/api/files/b/gb");
      else {
        expect(
          screen.queryByRole("link", { name: "下载 GenBank" }),
        ).not.toBeInTheDocument();
        expect(screen.getByText(/产物校验失败/)).toBeInTheDocument();
      }
    },
  );
});
