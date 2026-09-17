import { segmentInstance } from "./componentOrder";
import { palette } from "./palette";
import type { Feature, Preview } from "./types";

export type AnnotationRole = "promoter" | "rbs" | "gene" | "terminator";
export type AnnotationItem = { key: string; feature: Feature; role: AnnotationRole; name: string; group: string | null };
export type CassetteGroup = { key: string; feature: Feature; name: string; items: AnnotationItem[] };
export type AnnotationModel = { items: AnnotationItem[]; groups: CassetteGroup[]; remaining: Feature[] };
export type ReservedLabel = { x: number; y: number };

export const annotationRoles: AnnotationRole[] = ["promoter", "rbs", "gene", "terminator"];
export const annotationNames: Record<AnnotationRole, string> = {
  promoter: "启动子", rbs: "RBS", gene: "CDS", terminator: "终止子",
};
export const annotationColors: Record<AnnotationRole, string> = {
  promoter: "#d29922", rbs: "#ff7b72", gene: palette.gene, terminator: palette.terminator,
};
export const annotationTracks: Record<AnnotationRole, number> = {
  promoter: 116, rbs: 88, gene: 103, terminator: 94,
};

export function annotationRole(feature: Feature): AnnotationRole | null {
  const type = feature.type.toLowerCase();
  if (type === "promoter") return "promoter";
  if (type === "rbs") return "rbs";
  if (type === "terminator") return "terminator";
  if (type === "cds" || type === "gene") return "gene";
  return null;
}

export function annotationName(feature: Feature): string {
  const role = annotationRole(feature);
  if (role !== "gene") return role ? annotationNames[role] : feature.label;
  return feature.label.replace(/^cassette_\d+_/i, "").replace(/_CDS$/i, "") || feature.label;
}

export function annotationTooltip(feature: Feature): string {
  const strand = feature.strand === 1 ? "正链 →" : feature.strand === -1 ? "反链 ←" : "方向未标注";
  return `${feature.label} · ${feature.type} · ${feature.start_bp}–${feature.end_bp} bp · ${feature.end_bp - feature.start_bp + 1} bp · ${strand}`;
}

export function buildAnnotationModel(features: Feature[], segments: Preview["segments"], length: number): AnnotationModel {
  const valid = (f: Feature) => Number.isInteger(f.start_bp) && Number.isInteger(f.end_bp) && f.start_bp >= 1 && f.end_bp >= f.start_bp && f.end_bp <= length;
  const contains = (outer: { start_bp: number; end_bp: number }, f: Feature) => f.start_bp >= outer.start_bp && f.end_bp <= outer.end_bp;
  const expressions = segments.filter(s => s.kind === "expression");
  const owner = (f: Feature) => expressions.find(s => contains(s, f) && (!f.instance_id || f.instance_id === segmentInstance(s)));
  const groups: CassetteGroup[] = [];
  features.forEach((f, index) => {
    const match = /^expression_cassette_(\d+)$/i.exec(f.label);
    const segment = owner(f);
    if (match && f.type === "misc_feature" && valid(f) && segment) {
      groups.push({ key: `${segmentInstance(segment)}:${index}`, feature: f, name: `表达盒${match[1]}`, items: [] });
    }
  });
  const items: AnnotationItem[] = [];
  features.forEach((f, index) => {
    const role = annotationRole(f);
    if (!role || !valid(f) || (role !== "gene" && !owner(f))) return;
    // Gene and CDS commonly describe the same locus in imported GenBank files.
    if (f.type.toLowerCase() === "gene" && features.some(other => other.type.toLowerCase() === "cds" && other.start_bp === f.start_bp && other.end_bp === f.end_bp && other.strand === f.strand)) return;
    const group = groups.find(g => contains(g.feature, f));
    const item = { key: `feature-${index}`, feature: f, role, name: annotationName(f), group: group?.key || null };
    items.push(item);
    group?.items.push(item);
  });
  const grouped = new Set(groups.flatMap(g => [g.feature, ...g.items.map(i => i.feature)]));
  return { items, groups, remaining: features.filter(f => !grouped.has(f)) };
}

export function annotationAngle(bp: number, length: number) {
  return -Math.PI / 2 + ((bp - 1) / Math.max(1, length)) * Math.PI * 2;
}
export function annotationPoint(radius: number, angle: number) {
  return [220 + radius * Math.cos(angle), 220 + radius * Math.sin(angle)];
}
function pointString(radius: number, angle: number) {
  return annotationPoint(radius, angle).join(" ");
}
function curve(radius: number, start: number, end: number) {
  // Split long arcs, including full circles, to avoid coincident SVG endpoints.
  const pieces = Math.max(1, Math.ceil(Math.abs(end - start) / Math.PI));
  return Array.from({ length: pieces }, (_, i) =>
    `A ${radius} ${radius} 0 0 ${end >= start ? 1 : 0} ${pointString(radius, start + (end - start) * (i + 1) / pieces)}`,
  ).join(" ");
}
export function annotationArc(radius: number, start: number, end: number) {
  return `M ${pointString(radius, start)} ${curve(radius, start, end)}`;
}
export function genePath(feature: Feature, length: number) {
  const start = annotationAngle(feature.start_bp, length), end = annotationAngle(feature.end_bp + 1, length);
  const head = Math.min((end - start) / 3, 8 / 103);
  if (feature.strand === 1) {
    return `M ${pointString(109, start)} ${curve(109, start, end - head)} L ${pointString(103, end)} L ${pointString(97, end - head)} ${curve(97, end - head, start)} Z`;
  }
  if (feature.strand === -1) {
    return `M ${pointString(103, start)} L ${pointString(109, start + head)} ${curve(109, start + head, end)} L ${pointString(97, end)} ${curve(97, end, start + head)} Z`;
  }
  return `M ${pointString(109, start)} ${curve(109, start, end)} L ${pointString(97, end)} ${curve(97, end, start)} Z`;
}

export function textArc(radius: number, start: number, end: number) {
  return Math.sin((start + end) / 2) > 0 ? annotationArc(radius, end, start) : annotationArc(radius, start, end);
}
export function textWidth(text: string) {
  return Array.from(text).reduce((width, char) => width + (char.charCodeAt(0) > 127 ? 9 : 5.6), 0);
}

export function geneLabels(items: AnnotationItem[], length: number, reserved: ReservedLabel[]) {
  const occupied = reserved.map(label => ({ ...label }));
  const genes = items.filter(item => item.role === "gene");
  const labels = genes.map(item => {
    const start = annotationAngle(item.feature.start_bp, length), end = annotationAngle(item.feature.end_bp + 1, length);
    const mid = (start + end) / 2;
    const fits = (end - start) * 103 >= textWidth(item.name) + 20;
    if (fits) {
      const [x, y] = annotationPoint(103, mid);
      occupied.push({ x, y });
    }
    return { item, start, end, mid, fits };
  });
  return labels.map(label => {
    if (label.fits) return { ...label, leader: null };
    const side = Math.cos(label.mid) >= 0 ? 1 : -1;
    const x = side > 0 ? 412 : 28;
    const preferred = Math.max(44, Math.min(396, 220 + 176 * Math.sin(label.mid)));
    const candidates = Array.from({ length: 16 }, (_, i) => [preferred + i * 24, preferred - i * 24]).flat();
    const y = candidates.find(y => y >= 44 && y <= 396 && occupied.every(p => Math.abs(p.x - x) > 160 || Math.abs(p.y - y) >= 24));
    if (y === undefined) return { ...label, leader: null }; // Keep dense maps readable; full names remain in tooltips/details.
    occupied.push({ x, y });
    let name = label.item.name;
    while (textWidth(name) > 105 && name.length > 1) name = name.slice(0, -1);
    if (name !== label.item.name) name += "…";
    return { ...label, leader: { x, y, side, name } };
  });
}
