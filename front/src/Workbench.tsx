import { PlasmidRing } from "./PlasmidRing";
export { PlasmidRing } from "./PlasmidRing";
import {
  DragEvent,
  useCallback,
  useEffect,
  useMemo,
  useRef,
  useState,
} from "react";
import type { MouseEvent as ReactMouseEvent } from "react";

import type {
  Module,
  Context,
  Issue,
  Preview,
  Job,
  ComponentOrder,
  ComponentType,
} from "./types";
import {
  DEFAULT_ORDER,
  insertComponent,
  normalizeOrder,
  reorderGeometry,
  sameOrder,
  selectionGeometry,
} from "./componentOrder";
import { palette } from "./palette";
import { api } from "./api";
import { moduleName } from "./display";
type ModuleRole = Exclude<ComponentType, "expression">;
const moduleRole = (item: Module): ModuleRole =>
  item.type === "terminator" ? item.role! : item.type;
const slotTitle = (type: ModuleRole) =>
  type === "resistance"
    ? "抗性标记"
    : type === "replication"
      ? "复制起始位点"
      : type.toUpperCase();
function IssueList({ issues }: { issues: Issue[] }) {
  return (
    <>
      {issues.map((issue, index) => (
        <div className="issue" key={`${issue.code}-${index}`}>
          {issue.message}
        </div>
      ))}
    </>
  );
}

export function Workbench() {
  const [library, setLibrary] = useState<{
    resistance: Module[];
    replication: Module[];
    terminator: Module[];
  }>({ resistance: [], replication: [], terminator: [] });
  const [context, setContext] = useState<Context | null>(null);
  const [chosen, setChosen] = useState<{
    resistance_id: string;
    replication_id: string;
    t0_id: string | null;
    t1_id: string | null;
  }>({ resistance_id: "", replication_id: "", t0_id: null, t1_id: null });
  const [componentOrder, setComponentOrder] = useState<ComponentOrder>([
    ...DEFAULT_ORDER,
  ]);
  const [preview, setPreview] = useState<Preview | null>(null);
  const [geometry, setGeometry] = useState<Preview | null>(null);
  const [message, setMessage] = useState("正在读取项目…");
  const [job, setJob] = useState<Job | null>(null);
  const [collapsed, setCollapsed] = useState<Record<string, boolean>>({});
  const [filter, setFilter] = useState("");
  const [dragged, setDragged] = useState<Module | null>(null);
  const [dragTarget, setDragTarget] = useState<string | null>(null);
  const requestVersion = useRef(0);
  const loadVersion = useRef(0);
  const contextRef = useRef<Context | null>(null);
  const submitting = useRef(false);
  const ignoredJobs = useRef(new Set<string>());
  const [listDrag, setListDrag] = useState<{
    type: ComponentType;
    insertion: number;
    moved: boolean;
  } | null>(null);
  const listDragRef = useRef(listDrag);
  const listOrigin = useRef({ x: 0, y: 0 });
  const listRef = useRef<HTMLDivElement>(null);
  const [confirmWarnings, setConfirmWarnings] = useState<string[] | null>(null);
  const dialogRef = useRef<HTMLDivElement>(null);
  const generateButtonRef = useRef<HTMLButtonElement>(null);
  const [pending, setPending] = useState(false);
  const [refreshKey, setRefreshKey] = useState(0);
  const selected = useMemo(
    () => ({
      resistance: library.resistance.find(
        (item) => item.id === chosen.resistance_id,
      ),
      t0: library.terminator.find(
        (item) => item.role === "t0" && item.id === chosen.t0_id,
      ),
      t1: library.terminator.find(
        (item) => item.role === "t1" && item.id === chosen.t1_id,
      ),
      replication: library.replication.find(
        (item) => item.id === chosen.replication_id,
      ),
    }),
    [library, chosen],
  );
  const invalidateJob = useCallback(() => {
    setConfirmWarnings(null);
    if (contextRef.current?.active_job)
      ignoredJobs.current.add(contextRef.current.active_job.id);
    setJob((current) => {
      if (current) ignoredJobs.current.add(current.id);
      return null;
    });
  }, []);
  const choose = useCallback(
    (item: Module) => {
      requestVersion.current += 1;
      setChosen((current) => ({
        ...current,
        [`${moduleRole(item)}_id`]: item.id,
      }));
      setPreview(null);
      setGeometry(null);
      invalidateJob();
      setMessage("选择已更新，正在验证…");
    },
    [invalidateJob],
  );
  const reorder = useCallback(
    (nextOrder: ComponentOrder) => {
      if (sameOrder(nextOrder, componentOrder)) return;
      requestVersion.current += 1;
      setComponentOrder(nextOrder);
      setPreview(null);
      invalidateJob();
      setMessage("顺序已更新，等待验证…");
    },
    [componentOrder, invalidateJob],
  );
  const load = useCallback(async (preserve = false) => {
    const token = ++loadVersion.current;
    try {
      const [modules, nextContext] = await Promise.all([
        api<{
          resistance: Module[];
          replication: Module[];
          terminator?: Module[];
        }>("/api/modules"),
        api<Context>("/api/context"),
      ]);
      if (token !== loadVersion.current) return;
      const previous = contextRef.current;
      const changed =
        !previous ||
        previous.ready !== nextContext.ready ||
        previous.project.manifest_revision !==
          nextContext.project.manifest_revision ||
        previous.project.source_fingerprint !==
          nextContext.project.source_fingerprint;
      setLibrary({ ...modules, terminator: modules.terminator || [] });
      if (changed) {
        setConfirmWarnings(null);
        requestVersion.current += 1;
        setPreview(null);
        setGeometry(null);
      }
      contextRef.current = nextContext;
      setContext((current) =>
        changed ? nextContext : { ...nextContext, project: current!.project },
      );
      if (!preserve) {
        let saved: {
          resistance_id?: string;
          replication_id?: string;
          t0_id?: string | null;
          t1_id?: string | null;
          component_order?: ComponentOrder;
        } | null = null;
        try {
          saved = JSON.parse(
            localStorage.getItem(`plasmid:${nextContext.project.target}`) ||
              "null",
          );
        } catch {
          /* Storage may be disabled. */
        }
        const candidate = saved || nextContext.selection || {};
        setComponentOrder(normalizeOrder(candidate.component_order));
        setChosen({
          t0_id: (modules.terminator || []).some(
            (item) => item.role === "t0" && item.id === candidate.t0_id,
          )
            ? candidate.t0_id!
            : null,
          t1_id: (modules.terminator || []).some(
            (item) => item.role === "t1" && item.id === candidate.t1_id,
          )
            ? candidate.t1_id!
            : null,
          resistance_id: modules.resistance.some(
            (item) => item.id === candidate.resistance_id,
          )
            ? candidate.resistance_id || ""
            : "",
          replication_id: modules.replication.some(
            (item) => item.id === candidate.replication_id,
          )
            ? candidate.replication_id || ""
            : "",
        });
        setJob(nextContext.active_job);
      } else {
        setJob((current) => {
          if (
            nextContext.active_job &&
            !ignoredJobs.current.has(nextContext.active_job.id)
          )
            return nextContext.active_job;
          // A completed response is provisional until context reaches its commit.
          // Thereafter context owns both replacement artifacts and invalidation.
          if (
            current?.status === "succeeded" &&
            current.result &&
            nextContext.project.manifest_revision >=
              current.result.manifest_revision
          ) {
            return { ...current, result: nextContext.result || undefined };
          }
          return current;
        });
      }
      if (changed)
        setMessage(
          nextContext.ready
            ? "源上下文已读取，请选择模块并验证。"
            : "源构建尚未就绪；仍可查看可用模块和准备说明。",
        );
    } catch (error) {
      if (token === loadVersion.current)
        setMessage(`无法加载工作台：${(error as Error).message}`);
    }
  }, []);
  useEffect(() => {
    void load();
    const timer = window.setInterval(() => void load(true), 5000);
    return () => {
      window.clearInterval(timer);
      loadVersion.current += 1;
      requestVersion.current += 1;
    };
  }, [load]);
  useEffect(() => {
    if (context) {
      try {
        localStorage.setItem(
          `plasmid:${context.project.target}`,
          JSON.stringify({ ...chosen, component_order: componentOrder }),
        );
      } catch {
        /* Optional persistence. */
      }
    }
  }, [chosen, componentOrder, context]);
  const revision = context?.project.manifest_revision;
  const fingerprint = context?.project.source_fingerprint;
  const ready = context?.ready;
  useEffect(() => {
    if (!ready || !chosen.resistance_id || !chosen.replication_id) return;
    const token = ++requestVersion.current;
    void api<Preview>("/api/preview", {
      method: "POST",
      body: JSON.stringify({
        ...chosen,
        component_order: componentOrder,
        expected_revision: revision,
        source_fingerprint: fingerprint,
      }),
    })
      .then((next) => {
        if (token !== requestVersion.current) return;
        if (!sameOrder(next.component_order, componentOrder)) return;
        setPreview(next);
        setGeometry(next);
        setMessage(
          next.valid ? "预览已验证，可生成设计。" : "预览被验证规则阻止。",
        );
      })
      .catch((error) => {
        if (token !== requestVersion.current) return;
        setPreview(null);
        if (error.status === 409) {
          requestVersion.current += 1;
          setMessage("源构建已更新，正在重新读取上下文…");
          void load(true).then(() => setRefreshKey((value) => value + 1));
        } else setMessage(`预览失败：${error.message}`);
      });
    return () => {
      requestVersion.current += 1;
    };
  }, [chosen, componentOrder, revision, fingerprint, ready, load, refreshKey]);
  useEffect(() => {
    if (!job || !["queued", "running"].includes(job.status)) return;
    let cancelled = false;
    const timer = window.setTimeout(async () => {
      try {
        const next = await api<Job>(`/api/jobs/${job.id}`);
        if (cancelled) return;
        setJob(next);
        if (next.status === "succeeded") {
          setMessage("设计已完成，可以下载实际产物。");
          void load(true);
        }
        if (next.status === "failed")
          setMessage(`生成失败：${next.error?.message || next.message}`);
      } catch (error) {
        if (!cancelled) {
          setMessage(`任务状态读取失败：${(error as Error).message}`);
          setJob({ ...job });
        }
      }
    }, 1100);
    return () => {
      cancelled = true;
      window.clearTimeout(timer);
    };
  }, [job, load]);
  const generate = async (confirmed = false) => {
    if (
      !context?.ready ||
      !preview?.valid ||
      submitting.current ||
      (job && ["queued", "running"].includes(job.status))
    )
      return;
    if (!confirmed && (preview.terminator_warnings || []).length) {
      setConfirmWarnings(preview.terminator_warnings!);
      return;
    }
    setConfirmWarnings(null);
    submitting.current = true;
    setPending(true);
    const token = requestVersion.current;
    try {
      const next = await api<Job>("/api/generate", {
        method: "POST",
        body: JSON.stringify({
          ...chosen,
          component_order: componentOrder,
          expected_revision: revision,
          source_fingerprint: fingerprint,
        }),
      });
      if (token !== requestVersion.current) return;
      setJob(next);
      setMessage(`生成任务已提交：${next.stage}`);
      if (next.status === "succeeded") void load(true);
    } catch (error) {
      if (token !== requestVersion.current) return;
      setMessage(`生成失败：${(error as Error).message}`);
      if ((error as { status?: number }).status === 409) {
        requestVersion.current += 1;
        setPreview(null);
        void load(true).then(() => setRefreshKey((value) => value + 1));
      }
    } finally {
      submitting.current = false;
      setPending(false);
    }
  };
  const clear = (type: ModuleRole) => {
    requestVersion.current += 1;
    setPreview(null);
    setGeometry(null);
    invalidateJob();
    setChosen((current) => ({
      ...current,
      [`${type}_id`]: type === "t0" || type === "t1" ? null : "",
    }));
    setMessage("请选择两个模块以生成预览。");
  };
  const possibleResult = job?.result || context?.result;
  const currentResult =
    possibleResult &&
    context?.ready &&
    possibleResult.resistance_id === chosen.resistance_id &&
    possibleResult.replication_id === chosen.replication_id &&
    (possibleResult.t0_id ?? null) === chosen.t0_id &&
    (possibleResult.t1_id ?? null) === chosen.t1_id &&
    sameOrder(possibleResult.component_order, componentOrder) &&
    possibleResult.source_fingerprint === context.project.source_fingerprint
      ? possibleResult
      : null;
  const items = (type: Module["type"]) =>
    library[type].filter((item) =>
      `${item.name} ${item.antibiotic || ""} ${item.host_range || ""}`
        .toLowerCase()
        .includes(filter.toLowerCase()),
    );
  const endDrag = () => {
    setDragged(null);
    setDragTarget(null);
  };
  const startDrag = (event: DragEvent, item: Module) => {
    event.dataTransfer.setData("module-id", item.id);
    event.dataTransfer.setData("text/plain", item.id);
    event.dataTransfer.effectAllowed = "copy";
    setDragged(item);
  };
  const dragOver = (event: DragEvent, target: string, type?: ModuleRole) => {
    event.stopPropagation();
    if (type && dragged && moduleRole(dragged) !== type) {
      event.dataTransfer.dropEffect = "none";
      setDragTarget(null);
      return;
    }
    if (!dragged && !event.dataTransfer.types.includes("module-id")) return;
    event.preventDefault();
    event.dataTransfer.dropEffect = "copy";
    setDragTarget(target);
  };
  const dragLeave = (event: DragEvent) => {
    if (!event.currentTarget.contains(event.relatedTarget as Node | null))
      setDragTarget(null);
  };
  const drop = (event: DragEvent, type?: ModuleRole) => {
    event.preventDefault();
    event.stopPropagation();
    const id =
      event.dataTransfer.getData("module-id") ||
      event.dataTransfer.getData("text/plain");
    const candidates = type
      ? type === "t0" || type === "t1"
        ? library.terminator.filter((item) => item.role === type)
        : library[type]
      : [...library.resistance, ...library.replication, ...library.terminator];
    const item = candidates.find((candidate) => candidate.id === id);
    endDrag();
    if (item) choose(item);
  };
  const busy = pending || (!!job && ["queued", "running"].includes(job.status));
  const primaryFile =
    currentResult?.files.find((file) => file.id === "final_genbank") ||
    currentResult?.files.find((file) => file.filename.endsWith(".gb"));
  const allIssues = [...(context?.issues || []), ...(preview?.issues || [])];
  const warnings = [
    ...new Set([
      ...(preview?.warnings || []),
      ...(currentResult?.warnings || []),
    ]),
  ];
  const failed = /失败|无法/.test(message);
  const status = failed
    ? "操作失败"
    : busy
      ? "正在生成"
      : currentResult
        ? "已生成"
        : preview?.valid
          ? "检查通过"
          : allIssues.length
            ? "需要处理"
            : "待选择组件";
  const componentCount =
    Number(!!context?.construct) +
    Number(!!selected.resistance) +
    Number(!!selected.replication) +
    Number(!!selected.t0) +
    Number(!!selected.t1);
  const ringPreview = useMemo(
    () =>
      geometry
        ? reorderGeometry(geometry, componentOrder)
        : selectionGeometry(
            context?.construct || null,
            selected,
            componentOrder,
          ),
    [geometry, componentOrder, context?.construct, selected],
  );
  const draggingList = !!listDrag;
  const listIdentity = JSON.stringify([
    componentOrder,
    chosen,
    fingerprint,
    context?.construct?.sequence_sha256,
    ready,
  ]);
  useEffect(() => {
    if (listDragRef.current) {
      listDragRef.current = null;
      setListDrag(null);
    }
  }, [listIdentity]);
  useEffect(() => {
    if (!draggingList) return;
    const move = (event: MouseEvent) => {
      const current = listDragRef.current;
      if (!current) return;
      const remaining = componentOrder.filter((type) => type !== current.type);
      const insertion = remaining.filter((type) => {
        const box = listRef.current
          ?.querySelector(`[data-component-type='${type}']`)
          ?.getBoundingClientRect();
        return box && event.clientY > box.top + box.height / 2;
      }).length;
      const moved =
        current.moved ||
        Math.hypot(
          event.clientX - listOrigin.current.x,
          event.clientY - listOrigin.current.y,
        ) > 4;
      const next = { ...current, insertion, moved };
      listDragRef.current = next;
      setListDrag(next);
    };
    const cancel = () => {
      listDragRef.current = null;
      setListDrag(null);
    };
    const release = (event: MouseEvent) => {
      move(event);
      const current = listDragRef.current;
      cancel();
      if (current?.moved)
        reorder(
          insertComponent(componentOrder, current.type, current.insertion),
        );
    };
    const escape = (event: KeyboardEvent) => {
      if (event.key === "Escape") cancel();
    };
    window.addEventListener("mousemove", move);
    window.addEventListener("mouseup", release);
    window.addEventListener("keydown", escape);
    window.addEventListener("blur", cancel);
    return () => {
      window.removeEventListener("mousemove", move);
      window.removeEventListener("mouseup", release);
      window.removeEventListener("keydown", escape);
      window.removeEventListener("blur", cancel);
      listDragRef.current = null;
    };
  }, [draggingList]);
  const startListDrag = (event: ReactMouseEvent, type: ComponentType) => {
    if (
      event.button !== 0 ||
      (event.target as Element).closest("button") ||
      !(type === "expression" ? context?.construct : selected[type])
    )
      return;
    event.preventDefault();
    listOrigin.current = { x: event.clientX, y: event.clientY };
    const next = {
      type,
      insertion: componentOrder.indexOf(type),
      moved: false,
    };
    listDragRef.current = next;
    setListDrag(next);
  };
  const remainingList = componentOrder.filter(
    (type) => type !== listDrag?.type,
  );

  useEffect(() => {
    if (!confirmWarnings) return;
    dialogRef.current?.querySelector<HTMLButtonElement>("button")?.focus();
    const key = (event: KeyboardEvent) => {
      if (event.key === "Escape") {
        event.preventDefault();
        setConfirmWarnings(null);
        generateButtonRef.current?.focus();
      }
      if (event.key === "Tab") {
        const buttons = Array.from(
          dialogRef.current?.querySelectorAll<HTMLButtonElement>("button") ||
            [],
        );
        const first = buttons[0],
          last = buttons[buttons.length - 1];
        if (event.shiftKey && document.activeElement === first) {
          event.preventDefault();
          last?.focus();
        } else if (!event.shiftKey && document.activeElement === last) {
          event.preventDefault();
          first?.focus();
        }
      }
    };
    window.addEventListener("keydown", key);
    return () => window.removeEventListener("keydown", key);
  }, [confirmWarnings]);
  const cancelConfirmation = () => {
    setConfirmWarnings(null);
    generateButtonRef.current?.focus();
  };
  return (
    <div className="app">
      <header>
        <div className="header-title">
          <h1>🧬 质粒组装可视化</h1>
          <small>{context?.project.target || "正在连接"}</small>
        </div>
        <div className="header-actions">
          <button
            className="btn refresh"
            onClick={() => {
              requestVersion.current += 1;
              setPreview(null);
              setConfirmWarnings(null);
              void load(true).then(() => setRefreshKey((value) => value + 1));
            }}
          >
            刷新源构建
          </button>
          <button
            className="btn"
            onClick={() => {
              clear("resistance");
              clear("replication");
              clear("t0");
              clear("t1");
            }}
          >
            清空
          </button>
          <button
            ref={generateButtonRef}
            className={"btn generate" + (primaryFile ? "" : " btn-primary")}
            disabled={busy || !context?.ready || !preview?.valid}
            onClick={() => void generate()}
          >
            生成设计
          </button>
          {primaryFile && (
            <a
              className="btn btn-primary"
              aria-label={"下载 " + primaryFile.label}
              href={primaryFile.url}
              download={primaryFile.filename}
            >
              导出 .gb
            </a>
          )}
        </div>
      </header>
      <main>
        <aside className="library">
          <div className="pane-head">
            <h2>📦 组件库</h2>
          </div>
          <input
            className="search"
            aria-label="搜索模块"
            placeholder="搜索组件…"
            value={filter}
            onChange={(event) => setFilter(event.target.value)}
          />
          {(["replication", "resistance", "terminator"] as const).map(
            (type) => (
              <section className="category" key={type}>
                <button
                  className="category-title"
                  aria-expanded={!collapsed[type]}
                  onClick={() =>
                    setCollapsed((state) => ({
                      ...state,
                      [type]: !state[type],
                    }))
                  }
                >
                  <span>
                    <i style={{ background: palette[type] }} />
                    {type === "resistance"
                      ? "抗性模块"
                      : type === "replication"
                        ? "复制模块（ori）"
                        : "终止子"}
                  </span>
                  <span>{collapsed[type] ? "›" : "⌄"}</span>
                </button>
                {!collapsed[type] && (
                  <div className="category-items">
                    {items(type).length === 0 && (
                      <p className="empty">无匹配组件</p>
                    )}
                    {items(type).map((item) => (
                      <button
                        key={item.id}
                        className={
                          "module-card" +
                          (chosen[
                            (moduleRole(item) + "_id") as keyof typeof chosen
                          ] === item.id
                            ? " selected"
                            : "")
                        }
                        title={item.name}
                        draggable
                        onDragStart={(event) => startDrag(event, item)}
                        onDragEnd={endDrag}
                        onClick={() => choose(item)}
                        aria-label={"选择 " + item.name}
                      >
                        <i style={{ background: palette[type] }} />
                        <span className="module-name">{moduleName(item)}</span>
                        <small>{item.length_bp.toLocaleString()} bp</small>
                      </button>
                    ))}
                  </div>
                )}
              </section>
            ),
          )}
          <section className="category">
            <button
              className="category-title"
              aria-expanded={!collapsed.source}
              onClick={() =>
                setCollapsed((state) => ({ ...state, source: !state.source }))
              }
            >
              <span>
                <i style={{ background: palette.expression }} />
                表达构建
              </span>
              <span>{collapsed.source ? "›" : "⌄"}</span>
            </button>
            {!collapsed.source && (
              <div className="source-card" title={context?.construct?.name}>
                <i style={{ background: palette.expression }} />
                <span>
                  {context?.construct ? "完整表达构建" : "等待上游构建"}
                </span>
                <small>
                  {context?.construct?.length_bp.toLocaleString() || "—"} bp
                </small>
              </div>
            )}
          </section>
          <p className="note">点击组件，或拖到中间圆环、右侧列表。</p>
        </aside>
        <section
          className={"canvas" + (dragTarget === "canvas" ? " drop-target" : "")}
          onDragEnter={(event) => dragOver(event, "canvas")}
          onDragOver={(event) => dragOver(event, "canvas")}
          onDragLeave={dragLeave}
          onDrop={(event) => drop(event)}
        >
          {dragged && dragTarget === "canvas" && (
            <p className="drop-hint">
              松开以{selected[moduleRole(dragged)] ? "替换" : "添加"}{" "}
              {moduleName(dragged)}
            </p>
          )}
          <div className="construct-name">
            <h2>
              {context?.project.target
                ? context.project.target + " · 质粒预览"
                : "未命名质粒"}
            </h2>
            <div className="stats">
              <span>
                {preview
                  ? preview.length_bp.toLocaleString() + " bp"
                  : context?.construct
                    ? context.construct.length_bp.toLocaleString() + " bp"
                    : "0 bp"}
              </span>
              <span>组件：{componentCount}</span>
              <span>
                GC：
                {(
                  preview?.gc_percent ??
                  context?.construct?.gc_percent ??
                  0
                ).toFixed(2)}
                %
              </span>
            </div>
          </div>
          <PlasmidRing
            preview={ringPreview}
            construct={context?.construct || null}
            enzymes={context?.restriction_enzymes || []}
            componentOrder={componentOrder}
            onReorder={reorder}
          />
          {!preview && ringPreview && (
            <p className="geometry-pending">等待验证 · 当前组件示意</p>
          )}
          <div className="legend">
            {[
              ["replication", "复制模块"],
              ["resistance", "抗性模块"],
              ["expression", "表达构建"],
              ["terminator", "T0 / T1"],
            ].map(([kind, label]) => (
              <span key={kind}>
                <i style={{ background: palette[kind] }} />
                {label}
              </span>
            ))}
          </div>
        </section>
        <aside
          className={
            "assembly" + (dragTarget === "assembly" ? " drop-target" : "")
          }
          onDragEnter={(event) => dragOver(event, "assembly")}
          onDragOver={(event) => dragOver(event, "assembly")}
          onDragLeave={dragLeave}
          onDrop={(event) => drop(event)}
        >
          <div className="pane-head">
            <h2>📋 组装列表</h2>
            <span className="count">{componentCount}</span>
          </div>
          <div className="component-list" ref={listRef}>
            {componentOrder.map((type) => {
              const before =
                listDrag?.moved && remainingList[listDrag.insertion] === type;
              const after =
                listDrag?.moved &&
                listDrag.insertion === remainingList.length &&
                remainingList[remainingList.length - 1] === type;
              return (
                <div
                  key={type}
                  data-component-type={type}
                  className={
                    "component-row" +
                    (listDrag?.type === type && listDrag.moved
                      ? " sorting"
                      : "") +
                    (before ? " insertion-before" : "") +
                    (after ? " insertion-after" : "")
                  }
                  onMouseDown={(event) => startListDrag(event, type)}
                >
                  {type === "expression" ? (
                    <div
                      className="slot source-slot"
                      data-testid="source-slot"
                      title={context?.construct?.name}
                    >
                      <i style={{ background: palette.expression }} />
                      <b>
                        {context?.construct ? "完整表达构建" : "等待上游构建"}
                      </b>
                      <small>
                        {context?.construct?.length_bp.toLocaleString() || "—"}{" "}
                        bp
                      </small>
                      <span className="lock" title="由当前项目提供；可拖动排序">
                        🔒
                      </span>
                    </div>
                  ) : (
                    <Slot
                      title={slotTitle(type)}
                      item={selected[type]}
                      type={type}
                      onDrop={drop}
                      active={dragTarget === type}
                      onDragOver={dragOver}
                      onDragLeave={dragLeave}
                      onClear={clear}
                    />
                  )}
                </div>
              );
            })}
          </div>
          <section className="validation">
            <div className="status-line">
              <h3>检查结果</h3>
              <span
                className={
                  "status-badge" + (preview?.valid && !failed ? " ok" : "")
                }
                title={message}
              >
                {status}
              </span>
            </div>
            <p className="enzyme-summary">
              {(preview?.enzymes || context?.restriction_enzymes || [])
                .map((enzyme) => enzyme.name)
                .join(" / ") || "尚未指定组装酶"}
            </p>
            <IssueList issues={allIssues} />
            {failed && (
              <p className="issue" role="status">
                {message}
              </p>
            )}
            {busy && (
              <p className="job" role="status">
                {job?.message || "正在提交设计…"}
              </p>
            )}
          </section>
          {warnings.length > 0 && (
            <details className="foldout warnings">
              <summary>
                提示 <span>{warnings.length}</span>
              </summary>
              {warnings.map((warning) => (
                <p className="warning" key={warning}>
                  {warning}
                </p>
              ))}
            </details>
          )}
          <details className="foldout detail">
            <summary>组件详情与 DNA</summary>
            {[
              selected.resistance,
              selected.replication,
              selected.t0,
              selected.t1,
            ]
              .filter(Boolean)
              .map((item) => (
                <section className="module-detail" key={item!.id}>
                  <b>{moduleName(item!)}</b>
                  <p>{item!.name}</p>
                  <p>
                    {item!.antibiotic || item!.host_range} ·{" "}
                    {item!.resistance_gene || item!.copy_number}
                  </p>
                  <p>
                    {item!.length_bp.toLocaleString()} bp · GC{" "}
                    {item!.gc_percent.toFixed(2)}%
                  </p>
                  <details>
                    <summary>查看 DNA 序列</summary>
                    <textarea
                      aria-label={item!.name + " 只读 DNA"}
                      readOnly
                      value={item!.sequence}
                    />
                  </details>
                  {item!.notes.map((note) => (
                    <small key={note}>{note}</small>
                  ))}
                </section>
              ))}
            {!Object.values(selected).some(Boolean) && (
              <p className="note">选择组件后可查看详情。</p>
            )}
            {context?.construct && (
              <p className="note">
                {context.construct.name} ·{" "}
                {context.construct.length_bp.toLocaleString()} bp
              </p>
            )}
          </details>
          {currentResult && (
            <details className="foldout downloads">
              <summary>
                其他导出文件{" "}
                <span>
                  {currentResult.files.length - Number(!!primaryFile)}
                </span>
              </summary>
              {currentResult.files
                .filter((file) => file.id !== primaryFile?.id)
                .map((file) => (
                  <a key={file.id} href={file.url} download={file.filename}>
                    下载 {file.label}
                  </a>
                ))}
            </details>
          )}
        </aside>
      </main>
      {confirmWarnings && (
        <div
          className="terminator-backdrop"
          data-testid="terminator-backdrop"
          onClick={(event) => {
            if (event.target === event.currentTarget) cancelConfirmation();
          }}
        >
          <div
            className="terminator-dialog"
            ref={dialogRef}
            role="alertdialog"
            aria-modal="true"
            aria-labelledby="terminator-dialog-title"
            aria-describedby="terminator-dialog-warnings"
          >
            <h2 id="terminator-dialog-title">终止子配置提示</h2>
            <ul id="terminator-dialog-warnings">
              {confirmWarnings.map((warning, index) => (
                <li key={index}>{warning}</li>
              ))}
            </ul>
            <div className="dialog-actions">
              <button className="btn" onClick={cancelConfirmation}>
                返回调整
              </button>
              <button
                className="btn btn-primary"
                disabled={busy}
                onClick={() => void generate(true)}
              >
                继续生成
              </button>
            </div>
          </div>
        </div>
      )}
    </div>
  );
}

function Slot({
  title,
  item,
  type,
  onDrop,
  active,
  onDragOver,
  onDragLeave,
  onClear,
}: {
  title: string;
  item?: Module;
  type: ModuleRole;
  onDrop: (event: DragEvent, type: ModuleRole) => void;
  active: boolean;
  onDragOver: (event: DragEvent, target: string, type?: ModuleRole) => void;
  onDragLeave: (event: DragEvent) => void;
  onClear: (type: ModuleRole) => void;
}) {
  return (
    <div
      className={
        "slot" + (item ? "" : " empty-slot") + (active ? " drop-target" : "")
      }
      data-testid={type + "-slot"}
      onDragEnter={(event) => onDragOver(event, type, type)}
      onDragOver={(event) => onDragOver(event, type, type)}
      onDragLeave={onDragLeave}
      onDrop={(event) => onDrop(event, type)}
    >
      <i style={{ background: palette[type] }} />
      <b title={item?.name}>
        {item
          ? moduleName(item)
          : type === "t0" || type === "t1"
            ? "未添加" + title
            : "添加" + title}
      </b>
      {item && (
        <>
          <small>{item.length_bp.toLocaleString()} bp</small>
          <button
            className="remove"
            aria-label={"清空" + title}
            onClick={() => onClear(type)}
          >
            ×
          </button>
        </>
      )}
    </div>
  );
}
