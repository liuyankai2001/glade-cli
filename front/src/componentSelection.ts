import type { ComponentInstance, ComponentType, Module } from "./types";
import { COMPONENT_TYPES, normalizeOrder } from "./componentOrder";

export type Library = {
  resistance: Module[];
  replication: Module[];
  terminator: Module[];
  gap: Module[];
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
  return (parts || []).find(
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
  let componentIndexes: Array<number | null>;
  let originalIds: Set<string>;
  if ("components" in selection) {
    if (!Array.isArray(selection.components)) return initialComponents();
    const ids = new Set<string>();
    components = [];
    componentIndexes = [];
    for (const [index, raw] of selection.components.entries()) {
      if (
        !raw ||
        typeof raw !== "object" ||
        typeof raw.instance_id !== "string" ||
        !/^[A-Za-z0-9_.:-]{1,120}$/.test(raw.instance_id) ||
        ids.has(raw.instance_id) ||
        !COMPONENT_TYPES.includes(raw.component_type)
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
      componentIndexes.push(index);
    }
    if (
      components.filter((c) => c.component_type === "expression").length !==
        1 ||
      components.filter((c) => c.component_type === "replication").length !== 1
    )
      return initialComponents();
    originalIds = ids;
  } else {
    components = [];
    componentIndexes = [];
    originalIds = new Set<string>();
    let selectedIndex = 0;
    for (const type of normalizeOrder(selection.component_order)) {
      const moduleId = selection[`${type}_id`];
      const selected =
        type === "expression" || typeof moduleId === "string";
      if (!selected) {
        if (type === "replication") {
          components.push({ instance_id: type, component_type: type });
          componentIndexes.push(null);
          originalIds.add(type);
        }
        continue;
      }
      components.push({
        instance_id: type,
        component_type: type as ComponentType,
        ...(typeof moduleId === "string" ? { module_id: moduleId } : {}),
      });
      componentIndexes.push(selectedIndex);
      originalIds.add(type);
      selectedIndex += 1;
    }
  }
  const restored = components.flatMap((component, index) => {
    const originalIndex = componentIndexes[index];
    if (
      component.component_type === "expression" ||
      findModule(library, component)
    )
      return [{ component, originalIndex }];
    return component.component_type === "replication"
      ? [
          {
            component: {
              instance_id: component.instance_id,
              component_type: component.component_type,
            },
            originalIndex,
          },
        ]
      : [];
  });
  const isLegacy = selection.assembly_schema_version !== 2;
  if (
    !isLegacy ||
    components.some((component) => component.component_type === "gap") ||
    !(library.gap || []).length
  )
    return restored.map(({ component }) => component);

  const ids = new Set(originalIds);
  const gapFor = (component: ComponentInstance) =>
    library.gap.find((item) =>
      item.aliases?.includes(
        component.component_type === "resistance"
          ? "resistance_to_replication"
          : "replication_to_t1",
      ),
    );
  return restored.flatMap(({ component, originalIndex }) => {
    const gap =
      (component.component_type === "resistance" ||
        component.component_type === "replication") &&
      component.module_id &&
      gapFor(component);
    if (!gap || originalIndex === null) return [component];
    const base = `legacy-gap-${originalIndex}-${component.component_type}`;
    let instance_id = base;
    for (let suffix = 1; ids.has(instance_id); suffix += 1)
      instance_id = `${base}-${suffix}`;
    ids.add(instance_id);
    return [
      component,
      { instance_id, component_type: "gap", module_id: gap.id },
    ];
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
