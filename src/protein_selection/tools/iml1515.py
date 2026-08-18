"""Bounded, read-only queries over the bundled iML1515 model."""

from collections import defaultdict
import json
from pathlib import Path
import re
from typing import Any

from langchain_core.tools import BaseTool, StructuredTool

from src.config.config import GEM_DIR


DEFAULT_IML1515_PATH = GEM_DIR / "iML1515.json"
IML1515_TOOL_NAMES = (
    "find_iml1515_reactions",
    "get_iml1515_gene",
    "get_iml1515_reaction_context",
)

_GENE_TOKEN_PATTERN = re.compile(r"\bb\d{4}\b", re.IGNORECASE)
_GENE_ANNOTATION_FIELDS = (
    "uniprot",
    "ncbigene",
    "refseq_locus_tag",
    "refseq_name",
    "refseq_synonym",
)


class Iml1515ModelError(RuntimeError):
    """Raised when the local iML1515 model cannot be loaded safely."""


class Iml1515Repository:
    """In-memory, read-only index for the bundled COBRA JSON model."""

    def __init__(self, model_path: str | Path = DEFAULT_IML1515_PATH) -> None:
        self.model_path = Path(model_path)
        model = self._load_model(self.model_path)

        self.model_id = str(model.get("id", "iML1515"))
        self.model_version = str(model.get("version", "unknown"))
        self.compartments = dict(model.get("compartments", {}))
        self.reactions: list[dict[str, Any]] = model["reactions"]
        self.genes: list[dict[str, Any]] = model["genes"]

        self._reaction_by_id = {
            str(reaction["id"]).casefold(): reaction
            for reaction in self.reactions
        }
        self._gene_by_id = {
            str(gene["id"]).casefold(): gene for gene in self.genes
        }
        self._reaction_indexes = self._build_reaction_indexes()
        self._gene_indexes = self._build_gene_indexes()
        self._gene_reaction_ids = self._build_gene_reaction_index()

    @staticmethod
    def _load_model(model_path: Path) -> dict[str, Any]:
        try:
            with model_path.open("r", encoding="utf-8") as handle:
                model = json.load(handle)
        except FileNotFoundError as exc:
            raise Iml1515ModelError(
                f"iML1515 model file not found: {model_path}"
            ) from exc
        except (OSError, json.JSONDecodeError) as exc:
            raise Iml1515ModelError(
                f"failed to load iML1515 model: {model_path}"
            ) from exc

        if not isinstance(model, dict):
            raise Iml1515ModelError("iML1515 model root must be an object")
        for field in ("reactions", "genes"):
            if not isinstance(model.get(field), list):
                raise Iml1515ModelError(
                    f"iML1515 model field '{field}' must be a list"
                )
            for index, record in enumerate(model[field]):
                if not isinstance(record, dict) or not record.get("id"):
                    raise Iml1515ModelError(
                        f"iML1515 {field}[{index}] must contain an id"
                    )
        return model

    def _build_reaction_indexes(
        self,
    ) -> dict[str, dict[str, set[str]]]:
        indexes: dict[str, dict[str, set[str]]] = {
            "bigg": defaultdict(set),
            "kegg": defaultdict(set),
            "ec": defaultdict(set),
            "rhea": defaultdict(set),
        }
        for reaction in self.reactions:
            reaction_id = str(reaction["id"])
            indexes["bigg"][reaction_id.casefold()].add(reaction_id)
            annotation = reaction.get("annotation", {})
            for value in self._as_list(annotation.get("bigg.reaction")):
                indexes["bigg"][value.casefold()].add(reaction_id)
            for value in self._as_list(annotation.get("kegg.reaction")):
                indexes["kegg"][value.upper()].add(reaction_id)
            for value in self._as_list(annotation.get("ec-code")):
                indexes["ec"][self._normalize_ec(value)].add(reaction_id)
            for value in self._as_list(annotation.get("rhea")):
                indexes["rhea"][self._normalize_rhea(value)].add(reaction_id)
        return indexes

    def _build_gene_indexes(self) -> dict[str, set[str]]:
        indexes: dict[str, set[str]] = defaultdict(set)
        for gene in self.genes:
            gene_id = str(gene["id"])
            indexes[gene_id.casefold()].add(gene_id)
            name = gene.get("name")
            if name:
                indexes[str(name).casefold()].add(gene_id)
            annotation = gene.get("annotation", {})
            for field in _GENE_ANNOTATION_FIELDS:
                for value in self._as_list(annotation.get(field)):
                    indexes[value.casefold()].add(gene_id)
        return indexes

    def _build_gene_reaction_index(self) -> dict[str, list[str]]:
        index: dict[str, list[str]] = defaultdict(list)
        for reaction in self.reactions:
            reaction_id = str(reaction["id"])
            rule = str(reaction.get("gene_reaction_rule", ""))
            for token in self.extract_gene_tokens(rule):
                index[token].append(reaction_id)
        return index

    @staticmethod
    def _as_list(value: Any) -> list[str]:
        if value is None:
            return []
        if isinstance(value, list):
            return [str(item) for item in value]
        return [str(value)]

    @staticmethod
    def _normalize_ec(value: str) -> str:
        normalized = value.strip().upper()
        for prefix in ("EC:", "EC "):
            if normalized.startswith(prefix):
                normalized = normalized[len(prefix) :]
        return normalized

    @staticmethod
    def _normalize_rhea(value: str) -> str:
        normalized = value.strip().upper()
        if normalized.startswith("RHEA:"):
            normalized = normalized[5:]
        return normalized

    @staticmethod
    def extract_gene_tokens(gene_reaction_rule: str) -> list[str]:
        """Return unique b-number tokens without evaluating the GPR."""

        return list(
            dict.fromkeys(
                match.group(0).lower()
                for match in _GENE_TOKEN_PATTERN.finditer(
                    gene_reaction_rule
                )
            )
        )

    def find_reactions(
        self,
        *,
        bigg_reaction_id: str | None = None,
        kegg_reaction_id: str | None = None,
        ec_number: str | None = None,
        rhea_id: str | None = None,
        max_results: int = 20,
    ) -> dict[str, Any]:
        """Find reactions matching every supplied exact identifier."""

        if not 1 <= max_results <= 100:
            raise ValueError("max_results must be between 1 and 100")

        filters: list[set[str]] = []
        if bigg_reaction_id:
            filters.append(
                self._reaction_indexes["bigg"].get(
                    bigg_reaction_id.strip().casefold(),
                    set(),
                )
            )
        if kegg_reaction_id:
            filters.append(
                self._reaction_indexes["kegg"].get(
                    kegg_reaction_id.strip().upper(),
                    set(),
                )
            )
        if ec_number:
            filters.append(
                self._reaction_indexes["ec"].get(
                    self._normalize_ec(ec_number),
                    set(),
                )
            )
        if rhea_id:
            filters.append(
                self._reaction_indexes["rhea"].get(
                    self._normalize_rhea(rhea_id),
                    set(),
                )
            )
        if not filters:
            raise ValueError("at least one reaction identifier is required")

        matching_ids = set.intersection(*filters)
        ordered_ids = sorted(matching_ids)
        selected_ids = ordered_ids[:max_results]
        return {
            "model_id": self.model_id,
            "model_version": self.model_version,
            "match_count": len(ordered_ids),
            "truncated": len(ordered_ids) > len(selected_ids),
            "matches": [
                self._reaction_summary(self._reaction_by_id[item.casefold()])
                for item in selected_ids
            ],
        }

    def get_gene(
        self,
        identifier: str,
        *,
        max_results: int = 20,
    ) -> dict[str, Any]:
        """Find genes by locus, name, UniProt, NCBI, or RefSeq identifier."""

        normalized = identifier.strip().casefold()
        if not normalized:
            raise ValueError("gene identifier cannot be empty")
        if not 1 <= max_results <= 100:
            raise ValueError("max_results must be between 1 and 100")

        gene_ids = sorted(self._gene_indexes.get(normalized, set()))
        selected_ids = gene_ids[:max_results]
        return {
            "model_id": self.model_id,
            "query": identifier,
            "match_count": len(gene_ids),
            "truncated": len(gene_ids) > len(selected_ids),
            "matches": [
                self._gene_summary(self._gene_by_id[item.casefold()])
                for item in selected_ids
            ],
        }

    def get_reaction_context(self, bigg_reaction_id: str) -> dict[str, Any]:
        """Return an exact reaction, its raw GPR, genes, and compartments."""

        normalized = bigg_reaction_id.strip().casefold()
        if not normalized:
            raise ValueError("BiGG reaction identifier cannot be empty")
        reaction = self._reaction_by_id.get(normalized)
        if reaction is None:
            return {
                "model_id": self.model_id,
                "query": bigg_reaction_id,
                "found": False,
                "reaction": None,
                "genes": [],
                "unresolved_gene_tokens": [],
            }

        rule = str(reaction.get("gene_reaction_rule", ""))
        gene_tokens = self.extract_gene_tokens(rule)
        resolved_genes = [
            self._gene_summary(self._gene_by_id[token])
            for token in gene_tokens
            if token in self._gene_by_id
        ]
        unresolved = [
            token for token in gene_tokens if token not in self._gene_by_id
        ]
        compartment_ids = sorted(
            {
                metabolite_id.rsplit("_", 1)[-1]
                for metabolite_id in reaction.get("metabolites", {})
                if "_" in metabolite_id
            }
        )
        return {
            "model_id": self.model_id,
            "query": bigg_reaction_id,
            "found": True,
            "reaction": self._reaction_summary(reaction),
            "genes": resolved_genes,
            "unresolved_gene_tokens": unresolved,
            "compartments": {
                item: self.compartments.get(item, "unknown")
                for item in compartment_ids
            },
        }

    def _reaction_summary(self, reaction: dict[str, Any]) -> dict[str, Any]:
        rule = str(reaction.get("gene_reaction_rule", ""))
        return {
            "id": reaction.get("id"),
            "name": reaction.get("name"),
            "metabolites": reaction.get("metabolites", {}),
            "lower_bound": reaction.get("lower_bound"),
            "upper_bound": reaction.get("upper_bound"),
            "gene_reaction_rule": rule,
            "gene_tokens": self.extract_gene_tokens(rule),
            "subsystem": reaction.get("subsystem"),
            "annotation": reaction.get("annotation", {}),
        }

    def _gene_summary(self, gene: dict[str, Any]) -> dict[str, Any]:
        gene_id = str(gene["id"]).casefold()
        return {
            "id": gene.get("id"),
            "name": gene.get("name"),
            "annotation": gene.get("annotation", {}),
            "reaction_ids": self._gene_reaction_ids.get(gene_id, []),
        }


def build_iml1515_tools(
    model_path: str | Path = DEFAULT_IML1515_PATH,
) -> list[BaseTool]:
    """Build the three bounded LangChain tools over one cached model."""

    repository = Iml1515Repository(model_path)

    def find_iml1515_reactions(
        bigg_reaction_id: str | None = None,
        kegg_reaction_id: str | None = None,
        ec_number: str | None = None,
        rhea_id: str | None = None,
        max_results: int = 20,
    ) -> dict[str, Any]:
        return repository.find_reactions(
            bigg_reaction_id=bigg_reaction_id,
            kegg_reaction_id=kegg_reaction_id,
            ec_number=ec_number,
            rhea_id=rhea_id,
            max_results=max_results,
        )

    def get_iml1515_gene(
        identifier: str,
        max_results: int = 20,
    ) -> dict[str, Any]:
        return repository.get_gene(identifier, max_results=max_results)

    def get_iml1515_reaction_context(
        bigg_reaction_id: str,
    ) -> dict[str, Any]:
        return repository.get_reaction_context(bigg_reaction_id)

    return [
        StructuredTool.from_function(
            func=find_iml1515_reactions,
            name=IML1515_TOOL_NAMES[0],
            description=(
                "Find exact iML1515 reactions by BiGG, KEGG reaction, EC, "
                "or Rhea identifier. Multiple identifiers are ANDed."
            ),
        ),
        StructuredTool.from_function(
            func=get_iml1515_gene,
            name=IML1515_TOOL_NAMES[1],
            description=(
                "Find an iML1515 gene by b-number, gene name, UniProt, "
                "NCBI Gene, or RefSeq identifier."
            ),
        ),
        StructuredTool.from_function(
            func=get_iml1515_reaction_context,
            name=IML1515_TOOL_NAMES[2],
            description=(
                "Get one exact BiGG reaction with its raw GPR, associated "
                "genes, metabolites, annotations, and compartments."
            ),
        ),
    ]
