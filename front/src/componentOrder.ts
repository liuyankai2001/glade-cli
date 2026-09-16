import type {
  ComponentOrder,
  ComponentType,
  Context,
  Module,
  Preview,
} from "./types";

export const DEFAULT_ORDER: ComponentOrder = [
  "resistance",
  "replication",
  "t1",
  "expression",
  "t0",
];
export function isComponentOrder(value: unknown): value is ComponentOrder {
  return (
    Array.isArray(value) &&
    value.length === 5 &&
    new Set(value).size === 5 &&
    value.every((type) => DEFAULT_ORDER.includes(type))
  );
}
export function normalizeOrder(value: unknown): ComponentOrder {
  if (isComponentOrder(value)) return [...value];
  if (
    Array.isArray(value) &&
    value.length === 3 &&
    new Set(value).size === 3 &&
    value.every((type) =>
      ["resistance", "replication", "expression"].includes(type),
    )
  ) {
    return value.flatMap((type) =>
      type === "expression" ? ["t1", "expression", "t0"] : [type],
    ) as ComponentOrder;
  }
  return [...DEFAULT_ORDER];
}
export function sameOrder(a: unknown, b: unknown) {
  return normalizeOrder(a).join(",") === normalizeOrder(b).join(",");
}
export function insertComponent(
  order: ComponentOrder,
  type: ComponentType,
  index: number,
): ComponentOrder {
  const remaining = order.filter((item) => item !== type);
  remaining.splice(index, 0, type);
  return remaining;
}
export function segmentComponent(segment: {
  id: string;
  kind: string;
  component_type?: ComponentType;
}): ComponentType | null {
  if (segment.component_type && DEFAULT_ORDER.includes(segment.component_type))
    return segment.component_type;
  if (segment.id === "basic_seva_t0") return "t0";
  if (segment.id === "basic_seva_t1") return "t1";
  if (DEFAULT_ORDER.includes(segment.kind as ComponentType))
    return segment.kind as ComponentType;
  if (segment.id === "module_interval") return "resistance";
  if (["t1_interval", "t1"].includes(segment.id)) return "replication";
  if (["left_site", "right_site"].includes(segment.id)) return "expression";
  return null;
}
export function componentBlocks(preview: Pick<Preview, "segments">) {
  return DEFAULT_ORDER.flatMap((type) => {
    const parts = preview.segments.filter(
      (segment) => segmentComponent(segment) === type,
    );
    return parts.length
      ? [
          {
            type,
            start_bp: Math.min(...parts.map((part) => part.start_bp)),
            end_bp: Math.max(...parts.map((part) => part.end_bp)),
          },
        ]
      : [];
  }).sort((a, b) => a.start_bp - b.start_bp);
}

// Geometry only: carry server-owned blocks and annotations together while validation
// is pending. No client sequence is assembled or sent back to the server.
export function reorderGeometry(
  preview: Preview,
  order: ComponentOrder,
): Preview {
  const blocks = componentBlocks(preview);
  if (preview.segments.some((segment) => !segmentComponent(segment)))
    return preview;
  const offsets = new Map<ComponentType, number>();
  let cursor = 1;
  order.forEach((type) => {
    const block = blocks.find((block) => block.type === type);
    if (block) {
      offsets.set(type, cursor - block.start_bp);
      cursor += block.end_bp - block.start_bp + 1;
    }
  });
  const remap = (bp: number) => {
    const block = blocks.find(
      (block) => bp >= block.start_bp && bp <= block.end_bp,
    );
    return block ? bp + (offsets.get(block.type) || 0) : bp;
  };
  return {
    ...preview,
    valid: false,
    sequence: "",
    component_order: [...order],
    segments: preview.segments
      .map((segment) => ({
        ...segment,
        start_bp: remap(segment.start_bp),
        end_bp: remap(segment.end_bp),
      }))
      .sort((a, b) => a.start_bp - b.start_bp),
    features: preview.features.map((feature) => ({
      ...feature,
      start_bp: remap(feature.start_bp),
      end_bp: remap(feature.end_bp),
    })),
  };
}

export function selectionGeometry(
  construct: Context["construct"],
  selected: {
    resistance?: Module;
    replication?: Module;
    t0?: Module;
    t1?: Module;
  },
  order: ComponentOrder,
): Preview | null {
  if (!construct && !Object.values(selected).some(Boolean)) return null;
  let cursor = 1;
  let sourceOffset = 0;
  const segments: Preview["segments"] = [];
  order.forEach((type) => {
    const part = type === "expression" ? construct : selected[type];
    if (!part) return;
    if (type === "expression") sourceOffset = cursor - 1;
    segments.push({
      id: String(part.id),
      label: part.name,
      kind: type === "t0" || type === "t1" ? "terminator" : type,
      component_type: type,
      start_bp: cursor,
      end_bp: cursor + part.length_bp - 1,
      length_bp: part.length_bp,
    });
    cursor += part.length_bp;
  });
  return {
    valid: false,
    issues: [],
    warnings: [],
    terminator_warnings: [],
    sequence: "",
    length_bp: cursor - 1,
    gc_percent: 0,
    sequence_sha256: "",
    manifest_revision: 0,
    source_fingerprint: "",
    component_order: [...order],
    segments,
    features: (construct?.features || []).map((feature) => ({
      ...feature,
      start_bp: feature.start_bp + sourceOffset,
      end_bp: feature.end_bp + sourceOffset,
    })),
    enzymes: [],
  };
}
