import type { ComponentInstance, ComponentType, Module } from "./types";
import { DEFAULT_ORDER, normalizeOrder } from "./componentOrder";

export type Library = {
  resistance: Module[];
  replication: Module[];
  terminator: Module[];
};
export const initialComponents = (): ComponentInstance[] => [
  { instance_id: "replication", component_type: "replication" },
  { instance_id: "expression", component_type: "expression" },
];
export const moduleRole = (
  item: Module,
): Exclude<ComponentType, "expression"> =>
  item.type === "terminator" ? item.role! : item.type;

export function findModule(
  library: Library,
  component: ComponentInstance,
): Module | undefined {
  if (component.component_type === "expression") return undefined;
  const parts =
    component.component_type === "t0" || component.component_type === "t1"
      ? library.terminator
      : library[component.component_type];
  return parts.find(
    (item) =>
      item.id === component.module_id &&
      moduleRole(item) === component.component_type,
  );
}

export function restoreComponents(
  value: unknown,
  library: Library,
): ComponentInstance[] {
  if (!value || typeof value !== "object") return initialComponents();
  const selection = value as Record<string, unknown>;
  let components: ComponentInstance[];
  if ("components" in selection) {
    if (!Array.isArray(selection.components)) return initialComponents();
    const ids = new Set<string>();
    components = [];
    for (const raw of selection.components) {
      if (
        !raw ||
        typeof raw !== "object" ||
        typeof raw.instance_id !== "string" ||
        !/^[A-Za-z0-9_.:-]{1,120}$/.test(raw.instance_id) ||
        ids.has(raw.instance_id) ||
        !DEFAULT_ORDER.includes(raw.component_type)
      )
        return initialComponents();
      ids.add(raw.instance_id);
      components.push({
        instance_id: raw.instance_id,
        component_type: raw.component_type,
        ...(typeof raw.module_id === "string" &&
        raw.component_type !== "expression"
          ? { module_id: raw.module_id }
          : {}),
      });
    }
    if (
      components.filter((c) => c.component_type === "expression").length !==
        1 ||
      components.filter((c) => c.component_type === "replication").length !== 1
    )
      return initialComponents();
  } else {
    components = normalizeOrder(selection.component_order).map((type) => ({
      instance_id: type,
      component_type: type as ComponentType,
      ...(typeof selection[`${type}_id`] === "string"
        ? { module_id: selection[`${type}_id`] as string }
        : {}),
    }));
  }
  return components.flatMap((component) => {
    if (
      component.component_type === "expression" ||
      findModule(library, component)
    )
      return [component];
    return component.component_type === "replication"
      ? [
          {
            instance_id: component.instance_id,
            component_type: component.component_type,
          },
        ]
      : [];
  });
}

export function sameComponents(
  a: ComponentInstance[] | undefined,
  b: ComponentInstance[],
) {
  return (
    !!a &&
    a.length === b.length &&
    a.every(
      (item, index) =>
        item.instance_id === b[index].instance_id &&
        item.component_type === b[index].component_type &&
        item.module_id === b[index].module_id,
    )
  );
}
