export type ComponentType = "resistance" | "replication" | "expression";
export type ComponentOrder = ComponentType[];
export type Feature = {
  label: string;
  type: string;
  start_bp: number;
  end_bp: number;
  strand: number;
  kind: string;
};
export type Module = {
  id: string;
  name: string;
  type: "resistance" | "replication";
  sequence: string;
  length_bp: number;
  gc_percent: number;
  antibiotic?: string;
  resistance_gene?: string;
  host_range?: string;
  copy_number?: string;
  notes: string[];
};
export type Enzyme = {
  name: string;
  site: string;
  role: "left" | "right";
  overhang: string;
};
export type Context = {
  project: {
    target: string;
    name: string;
    host: string;
    manifest_revision: number;
    source_fingerprint: string;
  };
  ready: boolean;
  issues: Issue[];
  construct: {
    id: number;
    name: string;
    length_bp: number;
    gc_percent: number;
    sequence_sha256: string;
    features: Feature[];
  } | null;
  restriction_enzymes: Enzyme[];
  selection: {
    resistance_id: string;
    replication_id: string;
    component_order?: ComponentOrder;
  } | null;
  result: Result | null;
  active_job: Job | null;
};
export type Issue = {
  code: string;
  message: string;
  enzyme?: string;
  module?: string;
  start_bp?: number;
  end_bp?: number;
};
export type Preview = {
  component_order?: ComponentOrder;
  valid: boolean;
  issues: Issue[];
  warnings: string[];
  sequence: string;
  length_bp: number;
  gc_percent: number;
  sequence_sha256: string;
  segments: {
    id: string;
    label: string;
    kind: string;
    start_bp: number;
    end_bp: number;
    length_bp: number;
  }[];
  features: Feature[];
  enzymes: Enzyme[];
  source_fingerprint: string;
  manifest_revision: number;
};
export type Result = {
  component_order?: ComponentOrder;
  id: string;
  resistance_id: string;
  replication_id: string;
  sequence_sha256: string;
  length_bp: number;
  gc_percent: number;
  manifest_revision: number;
  source_fingerprint: string;
  files: {
    id: string;
    label: string;
    filename: string;
    format: string;
    url: string;
  }[];
  warnings: string[];
};
export type Job = {
  id: string;
  status: "queued" | "running" | "succeeded" | "failed";
  stage: string;
  message: string;
  error?: { code: string; message: string };
  result?: Result;
};
