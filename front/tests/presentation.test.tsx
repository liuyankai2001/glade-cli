import { afterEach, describe, expect, it } from "vitest";
import { cleanup, fireEvent, render, screen } from "@testing-library/react";
import "@testing-library/jest-dom/vitest";
import { PlasmidRing } from "../src/PlasmidRing";
import type { Preview } from "../src/types";

afterEach(cleanup);

describe("reference-style plasmid overview", () => {
  it("keeps gene annotations optional while preserving module and site coordinates", () => {
    const preview: Preview = {
      valid: true,
      issues: [],
      warnings: [],
      sequence: "ACGT",
      length_bp: 1000,
      gc_percent: 50,
      sequence_sha256: "sequence",
      manifest_revision: 1,
      source_fingerprint: "source",
      segments: [
        {
          id: "basic_seva_ap",
          label: "BASIC SEVA Ampicillin module with native regulation",
          kind: "resistance",
          start_bp: 1,
          end_bp: 400,
          length_bp: 400,
        },
      ],
      features: [
        {
          label: "protein_a",
          type: "CDS",
          kind: "gene",
          start_bp: 450,
          end_bp: 600,
          strand: 1,
        },
        {
          label: "EcoRI_retained_site",
          type: "misc_feature",
          kind: "restriction",
          start_bp: 401,
          end_bp: 406,
          strand: 1,
        },
      ],
      enzymes: [
        { name: "EcoRI", site: "GAATTC", role: "left", overhang: "AATT" },
      ],
    };
    const { container } = render(
      <PlasmidRing
        preview={preview}
        construct={null}
        enzymes={preview.enzymes}
      />,
    );
    expect(container.querySelector("[data-testid=strand-arrow]")).toBeNull();
    expect(
      container.querySelector("[data-testid=restriction-tick]"),
    ).toHaveAttribute("data-start-bp", "401");
    expect(screen.getByText("AmpR")).toBeInTheDocument();
    expect(screen.queryByText("EcoRI_retained_site")).not.toBeInTheDocument();
    fireEvent.click(screen.getByRole("button", { name: "显示基因注释" }));
    expect(
      container.querySelector("[data-testid=strand-arrow]"),
    ).toHaveAttribute("data-strand", "1");
    expect(container.querySelector("path[data-start-bp='1']")).toHaveAttribute(
      "data-end-bp",
      "400",
    );
  });
});
