import { useEffect, useRef, useState } from "react";
import type { MouseEvent as ReactMouseEvent } from "react";
import type {
  Preview,
  Context,
  Enzyme,
  ComponentOrder,
  ComponentType,
} from "./types";
import {
  componentBlocks,
  DEFAULT_ORDER,
  insertComponent,
  segmentComponent,
  segmentInstance,
} from "./componentOrder";
import { palette } from "./palette";
import { segmentName } from "./display";

type Segment = {
  id: string;
  label: string;
  kind: string;
  component_type?: ComponentType;
  instance_id?: string;
  start_bp: number;
  end_bp: number;
};

function arc(radius: number, start: number, end: number, width: number) {
  const point = (angle: number, r: number) => [
    220 + r * Math.cos(angle),
    220 + r * Math.sin(angle),
  ];
  const [x1, y1] = point(start, radius + width / 2);
  const [x2, y2] = point(end, radius + width / 2);
  const [x3, y3] = point(end, radius - width / 2);
  const [x4, y4] = point(start, radius - width / 2);
  const large = end - start > Math.PI ? 1 : 0;
  return `M ${x1} ${y1} A ${radius + width / 2} ${radius + width / 2} 0 ${large} 1 ${x2} ${y2} L ${x3} ${y3} A ${radius - width / 2} ${radius - width / 2} 0 ${large} 0 ${x4} ${y4} Z`;
}

export function PlasmidRing({
  preview,
  construct,
  componentOrder = DEFAULT_ORDER,
  onReorder,
}: {
  preview: Preview | null;
  construct: Context["construct"];
  enzymes: Enzyme[];
  componentOrder?: ComponentOrder;
  onReorder?: (order: ComponentOrder) => void;
}) {
  const [view, setView] = useState({ zoom: 1, rotate: 0 });
  const [annotations, setAnnotations] = useState(false);
  const svgRef = useRef<SVGSVGElement>(null);
  const [drag, setDrag] = useState<{
    type: string;
    angle: number;
    insertion: number;
    marker: number;
    moved: boolean;
  } | null>(null);
  const dragRef = useRef(drag);
  const origin = useRef({ x: 0, y: 0 });
  const length = preview?.length_bp || construct?.length_bp || 0;
  const outer: Segment[] =
    preview?.segments ||
    (construct
      ? [
          {
            id: "source",
            label: construct.name,
            kind: "expression",
            start_bp: 1,
            end_bp: construct.length_bp,
          },
        ]
      : []);
  const features = preview?.features || construct?.features || [];
  const angle = (bp: number) =>
    -Math.PI / 2 + ((Math.max(1, bp) - 1) / Math.max(1, length)) * Math.PI * 2;
  const genes = features.filter(
    (feature) => feature.type === "CDS" || feature.kind === "gene",
  );
  const sites = features.filter((feature) => feature.kind === "restriction");
  const blocks = componentBlocks({ segments: outer as Preview["segments"] });
  const dragging = !!drag;
  const geometryIdentity = JSON.stringify([
    length,
    componentOrder,
    view,
    preview?.source_fingerprint || construct?.sequence_sha256,
    outer.map((part) => [
      part.id,
      part.kind,
      part.component_type,
      part.instance_id,
      part.start_bp,
      part.end_bp,
    ]),
  ]);
  useEffect(() => {
    // Polling may deliver new objects with identical coordinates. Cancel only
    // when the geometry/source/view actually changes during a gesture.
    if (dragRef.current) {
      dragRef.current = null;
      setDrag(null);
    }
  }, [geometryIdentity]);
  useEffect(() => {
    if (!dragging) return;
    const position = (event: MouseEvent) => {
      const box = svgRef.current?.getBoundingClientRect();
      if (!box || !box.width || !box.height || !dragRef.current) return;
      // SVG defaults to xMidYMid meet; undo letterboxing and the centered view transform.
      const scale = Math.min(box.width, box.height) / 440;
      const x =
        (event.clientX - box.left - (box.width - scale * 440) / 2) / scale;
      const y =
        (event.clientY - box.top - (box.height - scale * 440) / 2) / scale;
      const a =
        Math.atan2((y - 220) / view.zoom, (x - 220) / view.zoom) -
        (view.rotate * Math.PI) / 180;
      const fraction =
        ((((a + Math.PI / 2) % (Math.PI * 2)) + Math.PI * 2) % (Math.PI * 2)) /
        (Math.PI * 2);
      const bp = fraction * length + 1;
      const remaining = blocks.filter(
        (block) => block.type !== dragRef.current!.type,
      );
      const insertion = remaining.filter(
        (block) => (block.start_bp + block.end_bp) / 2 < bp,
      ).length;
      const marker = remaining[insertion]?.start_bp || length + 1;
      const moved =
        dragRef.current.moved ||
        Math.hypot(
          event.clientX - origin.current.x,
          event.clientY - origin.current.y,
        ) > 4;
      const next = { ...dragRef.current, angle: a, insertion, marker, moved };
      dragRef.current = next;
      setDrag(next);
    };
    const cancel = () => {
      dragRef.current = null;
      setDrag(null);
    };
    const release = (event: MouseEvent) => {
      position(event);
      const current = dragRef.current;
      cancel();
      if (current?.moved) {
        const remaining = blocks.filter((block) => block.type !== current.type);
        // Hidden unselected groups retain their relative order; insert adjacent to
        // the visible group at the chosen boundary.
        const nextType = remaining[current.insertion]?.type;
        const allRemaining = componentOrder.filter(
          (type) => type !== current.type,
        );
        const index = nextType
          ? allRemaining.indexOf(nextType)
          : allRemaining.length;
        onReorder?.(insertComponent(componentOrder, current.type, index));
      }
    };
    const escape = (event: KeyboardEvent) => {
      if (event.key === "Escape") cancel();
    };
    window.addEventListener("mousemove", position);
    window.addEventListener("mouseup", release);
    window.addEventListener("keydown", escape);
    window.addEventListener("blur", cancel);
    return () => {
      window.removeEventListener("mousemove", position);
      window.removeEventListener("mouseup", release);
      window.removeEventListener("keydown", escape);
      window.removeEventListener("blur", cancel);
      dragRef.current = null;
    };
  }, [dragging]);
  const startDrag = (event: ReactMouseEvent, part: Segment) => {
    const type = segmentInstance(part);
    if (event.button !== 0 || !type || !onReorder) return;
    event.preventDefault();
    event.stopPropagation();
    origin.current = { x: event.clientX, y: event.clientY };
    const block = blocks.find((block) => block.type === type)!;
    const next = {
      type,
      angle: angle(block.start_bp),
      insertion: componentOrder.indexOf(type),
      marker: block.start_bp,
      moved: false,
    };
    dragRef.current = next;
    setDrag(next);
  };

  // Keep short independent components named even with gene annotations hidden.
  // Place leaders in side columns and keep their labels away from enzyme labels.
  const callouts = outer
    .filter((part) =>
      ["t0", "t1", "gap"].includes(segmentComponent(part) || ""),
    )
    .map((part) => {
      const a = (angle(part.start_bp) + angle(part.end_bp + 1)) / 2;
      const side = Math.cos(a) >= 0 ? 1 : -1;
      return {
        part,
        a,
        side,
        x: side > 0 ? 412 : 28,
        y: Math.max(44, Math.min(396, 220 + 176 * Math.sin(a))),
      };
    })
    .sort((a, b) => a.y - b.y);
  const occupied = sites.map((site) => ({
    x: 220 + 187 * Math.cos(angle(site.start_bp)),
    y: 220 + 187 * Math.sin(angle(site.start_bp)),
  }));
  callouts.forEach((callout) => {
    const free = (y: number) =>
      occupied.every(
        (label) =>
          Math.abs(label.x - callout.x) > 76 || Math.abs(label.y - y) >= 24,
      );
    const candidates = Array.from({ length: 16 }, (_, i) => [
      callout.y + i * 24,
      callout.y - i * 24,
    ]).flat();
    callout.y =
      candidates.find((y) => y >= 44 && y <= 396 && free(y)) ?? callout.y;
    occupied.push({ x: callout.x, y: callout.y });
  });
  const draw = (part: Segment, radius: number, width: number, label = true) => {
    const start = angle(part.start_bp),
      end = angle(part.end_bp + 1);
    const independentTerminator = ["t0", "t1"].includes(
      segmentComponent(part) || "",
    );
    const middle = (start + end) / 2;
    const visible =
      label &&
      !independentTerminator &&
      !["linker", "restriction"].includes(part.kind) &&
      ((part.end_bp - part.start_bp + 1) / Math.max(1, length)) * 360 > 22;
    return (
      <g
        key={`${part.instance_id || part.id}:${part.id}:${part.start_bp}`}
        data-ring-component={
          radius === 138 ? segmentComponent(part) || undefined : undefined
        }
        data-ring-instance={
          radius === 138 ? segmentInstance(part) || undefined : undefined
        }
        className={
          radius === 138 && onReorder && segmentComponent(part)
            ? "ring-component"
            : undefined
        }
        onMouseDown={
          radius === 138 ? (event) => startDrag(event, part) : undefined
        }
      >
        {part.end_bp - part.start_bp + 1 >= length ? (
          <circle
            data-testid="full-segment"
            cx="220"
            cy="220"
            r={radius}
            fill="none"
            stroke={palette[part.kind] || "#77839a"}
            strokeWidth={width}
          />
        ) : (
          <path
            data-start-bp={part.start_bp}
            data-end-bp={part.end_bp}
            d={arc(
              radius,
              start,
              independentTerminator ? Math.max(end, start + 0.015) : end,
              width,
            )}
            fill={palette[part.kind] || "#77839a"}
            stroke="#0d1117"
            strokeWidth={independentTerminator ? "0.4" : "1.5"}
          />
        )}
        <title>
          {part.label} · {part.start_bp}–{part.end_bp} bp
        </title>
        {visible && (
          <text
            x={220 + 176 * Math.cos(middle)}
            y={220 + 176 * Math.sin(middle)}
            textAnchor="middle"
          >
            {segmentName(part)}
          </text>
        )}
      </g>
    );
  };

  return (
    <div className="ring-wrap">
      <svg
        ref={svgRef}
        className="ring"
        viewBox="0 0 440 440"
        aria-label="质粒环图"
        role="img"
      >
        <g
          transform={`rotate(${view.rotate} 220 220) scale(${view.zoom}) translate(${220 / view.zoom - 220} ${220 / view.zoom - 220})`}
        >
          <circle
            cx="220"
            cy="220"
            r="138"
            fill="none"
            stroke="#273142"
            strokeWidth="28"
          />
          {outer.map((part) => draw(part, 138, 28))}
          {callouts.map(({ part, a, side, x, y }) => (
            <g
              key={`callout-${part.instance_id || part.id}`}
              className={
                segmentComponent(part) === "gap"
                  ? "gap-leader"
                  : "terminator-leader"
              }
              onMouseDown={(event) => startDrag(event, part)}
            >
              <path
                d={`M ${220 + 153 * Math.cos(a)} ${220 + 153 * Math.sin(a)} L ${x - side * 40} ${y - 4} L ${x - side * 4} ${y - 4}`}
                fill="none"
                stroke={palette[segmentComponent(part) || part.kind] || "#77839a"}
                strokeWidth="1"
              />
              <text
                className={
                  segmentComponent(part) === "gap"
                    ? "gap-callout"
                    : "terminator-callout"
                }
                style={{
                  fill:
                    palette[segmentComponent(part) || part.kind] || "#77839a",
                }}
                x={x}
                y={y - 8}
                textAnchor={side > 0 ? "end" : "start"}
              >
                {segmentName(part)}
              </text>
            </g>
          ))}
          {drag?.moved &&
            (() => {
              const block = blocks.find((block) => block.type === drag.type)!;
              const kind = segmentComponent(
                outer.find((part) => segmentInstance(part) === drag.type)!,
              );
              const span =
                ((block.end_bp - block.start_bp + 1) / length) * Math.PI * 2;
              const markerAngle = angle(drag.marker);
              return (
                <g className="ring-drag-overlay" pointerEvents="none">
                  <path
                    data-testid="ring-drag-ghost"
                    d={arc(
                      155,
                      drag.angle - span / 2,
                      drag.angle + span / 2,
                      28,
                    )}
                    fill={palette[kind || ""]}
                    opacity=".65"
                  />
                  <line
                    data-testid="ring-insertion-marker"
                    x1={220 + 116 * Math.cos(markerAngle)}
                    y1={220 + 116 * Math.sin(markerAngle)}
                    x2={220 + 178 * Math.cos(markerAngle)}
                    y2={220 + 178 * Math.sin(markerAngle)}
                    style={{
                      stroke: palette[kind || ""] || "#77839a",
                      strokeWidth: 4,
                    }}
                  />
                </g>
              );
            })()}
          {annotations &&
            genes.map((feature, index) => {
              const a = angle(
                feature.strand < 0 ? feature.start_bp : feature.end_bp + 1,
              );
              return (
                <g key={`feature-${index}`}>
                  {draw({ ...feature, id: `feature-${index}` }, 103, 12, false)}
                  {feature.strand !== 0 && (
                    <polygon
                      data-testid="strand-arrow"
                      data-strand={feature.strand}
                      points="-4,-4 4,0 -4,4"
                      fill="#e7edf6"
                      transform={`translate(${220 + 103 * Math.cos(a)} ${220 + 103 * Math.sin(a)}) rotate(${(a * 180) / Math.PI + (feature.strand > 0 ? 90 : -90)})`}
                    />
                  )}
                </g>
              );
            })}
          {sites.map((site, index) => {
            const a = angle(site.start_bp);
            return (
              <g key={`site-${index}`}>
                <line
                  data-testid="restriction-tick"
                  data-start-bp={site.start_bp}
                  x1={220 + 122 * Math.cos(a)}
                  y1={220 + 122 * Math.sin(a)}
                  x2={220 + 160 * Math.cos(a)}
                  y2={220 + 160 * Math.sin(a)}
                />
                <text
                  className="site-label"
                  x={220 + 187 * Math.cos(a)}
                  y={220 + 187 * Math.sin(a)}
                  textAnchor="middle"
                >
                  {segmentName(site)}
                </text>
              </g>
            );
          })}
        </g>
        <circle cx="220" cy="220" r="78" fill="#0d1117" />
        <text x="220" y="214" textAnchor="middle" className="ring-number">
          {length ? length.toLocaleString() : "—"}
        </text>
        <text x="220" y="236" textAnchor="middle" className="ring-unit">
          bp
        </text>
      </svg>
      <div className="ring-controls">
        <button
          aria-label="缩小环图"
          onClick={() =>
            setView((value) => ({
              ...value,
              zoom: Math.max(0.8, value.zoom - 0.1),
            }))
          }
        >
          −
        </button>
        <button
          aria-label="旋转环图"
          onClick={() =>
            setView((value) => ({ ...value, rotate: value.rotate + 30 }))
          }
        >
          ↻
        </button>
        <button
          aria-label="放大环图"
          onClick={() =>
            setView((value) => ({
              ...value,
              zoom: Math.min(1.15, value.zoom + 0.1),
            }))
          }
        >
          +
        </button>
        <button
          className="annotation-toggle"
          aria-label={annotations ? "隐藏基因注释" : "显示基因注释"}
          aria-pressed={annotations}
          onClick={() => setAnnotations((value) => !value)}
        >
          基因注释
        </button>
      </div>
      {annotations && (
        <details className="annotation-details">
          <summary>注释详情 · {features.length} 项</summary>
          <ul className="feature-legend">
            {features.map((feature, index) => (
              <li key={index}>
                <i style={{ background: palette[feature.kind] || "#77839a" }} />{" "}
                {feature.label} · {feature.start_bp}–{feature.end_bp} bp{" "}
                {feature.strand === 1 ? "→" : feature.strand === -1 ? "←" : ""}
              </li>
            ))}
          </ul>
        </details>
      )}
    </div>
  );
}
