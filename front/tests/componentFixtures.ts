import type { ComponentInstance, Preview } from "../src/types";

// Adapt the API double to the submitted ordered instances, preserving literal
// fixture lengths and coordinates within each module block.
export function componentPreview(
  preview: Preview,
  init?: RequestInit,
): Preview {
  const components: ComponentInstance[] | undefined = JSON.parse(
    String(init?.body || "{}"),
  ).components;
  if (!components) return preview;
  let cursor = 1;
  const segments = components.flatMap((component) => {
    const parts = preview.segments.filter((segment) => {
      const role =
        segment.component_type ||
        (segment.id === "module_interval"
          ? "resistance"
          : ["t1", "t1_interval"].includes(segment.id)
            ? "replication"
            : ["left_site", "right_site"].includes(segment.id)
              ? "expression"
              : segment.kind);
      return role === component.component_type;
    });
    const origin = Math.min(...parts.map((p) => p.start_bp));
    const offset = cursor - origin;
    if (parts.length)
      cursor += Math.max(...parts.map((p) => p.end_bp)) - origin + 1;
    return parts.map((segment) => ({
      ...segment,
      instance_id: component.instance_id,
      component_type: component.component_type,
      start_bp: segment.start_bp + offset,
      end_bp: segment.end_bp + offset,
    }));
  });
  return {
    ...preview,
    components,
    component_order: components.map((c) => c.instance_id),
    segments,
  };
}
