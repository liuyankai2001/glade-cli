const names: Record<string, string> = {
  basic_seva_ap: "AmpR",
  basic_seva_km: "KanR",
  basic_seva_cm: "CmR",
  basic_seva_sm_sp: "Sm/SpR",
  basic_seva_tet_5a: "TetR",
  basic_seva_gm: "GmR",
  basic_seva_gm_11: "GmR (66.11)",
  basic_seva_rsf1010: "RSF1010",
  basic_seva_p15a: "p15A",
  basic_seva_psc101: "pSC101",
  basic_seva_pbr322_rop: "pBR322 (rop)",
  basic_seva_psc101_pkd46_ts: "pSC101 (温敏)",
};

export function moduleName(item: { id: string; name: string }): string {
  return names[item.id] || item.name;
}

export function segmentName(item: {
  id?: string;
  label: string;
  kind: string;
}): string {
  if (item.id && names[item.id]) return names[item.id];
  if (item.kind === "expression") return "表达构建";
  if (item.kind === "restriction")
    return item.label.replace(/_retained_site$/, "");
  return item.label;
}
