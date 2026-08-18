"""Database-evidence subagent backed by the ToolUniverse MCP server."""

from collections.abc import Sequence
from typing import Literal, Self

from deepagents.middleware.subagents import SubAgent
from langchain.agents.middleware.types import AgentMiddleware
from langchain_core.language_models import BaseChatModel
from langchain_core.tools import BaseTool
from pydantic import BaseModel, ConfigDict, Field, model_validator

from src.protein_selection.agents.dependency_evidence import CuratedDependencyAssertion
from src.protein_selection.integrations.tooluniverse import (
    TOOLUNIVERSE_DEFAULT_TOOL_NAMES,
    ToolUniverseConfig,
    load_tooluniverse_tools,
)


BIO_DATABASE_RESEARCHER_NAME = "bio_database_researcher"
BIO_DATABASE_RESEARCHER_TOOL_NAMES = TOOLUNIVERSE_DEFAULT_TOOL_NAMES

EvidenceDatabase = Literal[
    "UniProt",
    "EBI Proteins",
    "KEGG",
    "Rhea",
    "Complex Portal",
    "IntAct",
    "STRING",
    "InterPro",
]
EvidenceType = Literal[
    "curated_annotation",
    "reaction_mapping",
    "curated_complex",
    "physical_interaction",
    "functional_association",
    "domain_inference",
]
EvidenceDirection = Literal["supports", "contradicts", "context_only"]
DependencyNecessity = Literal[
    "required",
    "enhancing",
    "associated",
    "uncertain",
]
AssessmentStatus = Literal["supported", "not_supported", "uncertain"]
ReactionMatch = Literal["matched", "mismatched", "uncertain"]


class DatabaseEvidence(BaseModel):
    """One traceable claim derived from a database record."""

    model_config = ConfigDict(extra="forbid", str_strip_whitespace=True)

    evidence_id: str = Field(
        min_length=1,
        description="Unique identifier within this report",
    )
    database: EvidenceDatabase
    record_id: str = Field(min_length=1)
    source_url: str | None = None
    claim: str = Field(min_length=1)
    evidence_type: EvidenceType
    direction: EvidenceDirection
    related_proteins: list[str] = Field(default_factory=list)
    limitations: list[str] = Field(default_factory=list)
    raw_evidence_ids: list[str] = Field(default_factory=list)
    source_locator: str | None = None
    supporting_excerpt: str | None = None


class ProteinDependencyCandidate(BaseModel):
    """A native protein partner indicated by database evidence."""

    model_config = ConfigDict(extra="forbid", str_strip_whitespace=True)

    protein_name: str = Field(min_length=1)
    uniprot_id: str | None = None
    organism_name: str | None = None
    taxon_id: int | None = None
    role: str = Field(min_length=1)
    necessity: DependencyNecessity
    evidence_ids: list[str] = Field(min_length=1)
    assertion_ids: list[str] = Field(default_factory=list)


class SmallMoleculeCofactor(BaseModel):
    """A non-protein cofactor kept separate from protein dependencies."""

    model_config = ConfigDict(extra="forbid", str_strip_whitespace=True)

    name: str = Field(min_length=1)
    database_id: str | None = None
    role: str = Field(min_length=1)
    evidence_ids: list[str] = Field(min_length=1)


class BioDatabaseResearchResult(BaseModel):
    """Structured evidence report returned by the database researcher."""

    model_config = ConfigDict(extra="forbid", str_strip_whitespace=True)

    input_uniprot_id: str = Field(min_length=1)
    reaction_id: str = Field(min_length=1)
    source_organism: str | None = None
    source_taxon_id: int | None = None
    reaction_match: ReactionMatch
    auxiliary_requirement_signal: AssessmentStatus
    native_system_summary: str = Field(min_length=1)
    candidate_protein_dependencies: list[ProteinDependencyCandidate] = Field(
        default_factory=list
    )
    small_molecule_cofactors: list[SmallMoleculeCofactor] = Field(
        default_factory=list
    )
    evidence: list[DatabaseEvidence] = Field(default_factory=list)
    dependency_assertions: list[CuratedDependencyAssertion] = Field(
        default_factory=list
    )
    conflicting_evidence_ids: list[str] = Field(default_factory=list)
    unresolved_questions: list[str] = Field(default_factory=list)

    @model_validator(mode="after")
    def validate_evidence_references(self) -> Self:
        """Ensure every evidence reference resolves within the report."""

        evidence_ids = [item.evidence_id for item in self.evidence]
        if len(evidence_ids) != len(set(evidence_ids)):
            raise ValueError("evidence_id values must be unique")

        known_ids = set(evidence_ids)
        assertion_ids = [
            item.assertion_id for item in self.dependency_assertions
        ]
        if len(assertion_ids) != len(set(assertion_ids)):
            raise ValueError("dependency assertion IDs must be unique")
        known_assertion_ids = set(assertion_ids)
        referenced_ids = set(self.conflicting_evidence_ids)
        for candidate in self.candidate_protein_dependencies:
            referenced_ids.update(candidate.evidence_ids)
            unknown_assertions = sorted(
                set(candidate.assertion_ids) - known_assertion_ids
            )
            if unknown_assertions:
                raise ValueError(
                    "unknown dependency assertion references: "
                    + ", ".join(unknown_assertions)
                )
        for cofactor in self.small_molecule_cofactors:
            referenced_ids.update(cofactor.evidence_ids)
        for assertion in self.dependency_assertions:
            referenced_ids.update(assertion.evidence_ids)

        unknown_ids = sorted(referenced_ids - known_ids)
        if unknown_ids:
            unknown_text = ", ".join(unknown_ids)
            raise ValueError(f"unknown evidence references: {unknown_text}")
        return self


BIO_DATABASE_RESEARCHER_PROMPT = """你是辅助蛋白检索流程中的生物数据库证据研究员。

输入包含一个经过校验的 UniProt 蛋白 ID 和 KEGG 反应 ID。输入蛋白可能来自任意
物种并被异源导入大肠杆菌。你的任务仅限于重建该蛋白在来源物种中的原生催化
系统和蛋白依赖证据；不要寻找大肠杆菌替代蛋白，也不要替最终裁决智能体决定
最终补全列表。

必须遵循以下调查顺序：
1. 首先查询 UniProt，确认来源物种、NCBI taxonomy ID、功能注释和可用交叉引用。
2. 查询 KEGG 目标反应及其 EC 映射，再用 Rhea 核对反应方向和参与物。一个反应
   对应多个 EC 或多个蛋白家族时，必须保留歧义。
3. 从 UniProt、Complex Portal 和 IntAct 查找明确的复合物组成、稳定亚基、电子
   传递伙伴以及成熟或装配蛋白。
4. STRING 只能提供功能关联线索；除非有独立的直接或人工整理证据，否则不能
   将 STRING 关联写成必需或增效蛋白。InterPro 结构域同样只能作为推断证据。
5. 对蛋白名称或基因名称搜索，必须显式传入来源物种或 taxonomy ID。STRING 和
   Complex Portal 工具存在默认人类物种，绝对不能依赖默认值。若数据库不覆盖
   来源物种，记录限制，不得用人类结果冒充原生证据。
6. 严格区分蛋白依赖与金属离子、辅酶、血红素等小分子辅因子。小分子只能写入
   small_molecule_cofactors。
7. 每个事实都建立独立 evidence_id，并记录数据库、记录 ID、可用 URL、证据类型
   和局限。候选蛋白必须引用这些 evidence_id。不得编造 UniProt ID、复合物成员、
   数据库字段、URL 或实验结论。

auxiliary_requirement_signal 的含义：
- supported：数据库存在直接或人工整理证据支持 required 或 enhancing
  蛋白伙伴；候选的 necessity 必须明确区分两者；
- not_supported：存在正面证据表明功能催化单元不包含另一种蛋白；
- uncertain：只有关联/结构域线索、数据库无覆盖、证据不足或互相冲突。

required 只能用于正面证据表明缺少伙伴时不能催化的情形。如果输入蛋白
单独仍有可检测催化活性，而伙伴仅提高活性、稳定性或调控，候选必须标为
enhancing。异源复合物组成和 iML1515 GPR 不能单独区分这两种情形。

若数据库记录逐字、明确地在完整蛋白层级断言某候选对指定反应 required 或
enhancing，填写 CuratedDependencyAssertion，并记录物种、反应 ID、独立 lineage
及逐字 EvidenceSpan。Complex Portal 组成、IntAct/STRING 互作、iML1515 GPR、
同一论文的数据库转引只能标为 composition/gpr/interaction，不能填写为
explicit_dependency，也不能单独把 auxiliary_requirement_signal 提升为 supported。

数据库没有记录不能作为“不需要辅助蛋白”的证据。最终输出必须符合指定的
BioDatabaseResearchResult 结构。只保留与候选依赖直接相关的最多 5 条 evidence；
无法填写完整 CuratedDependencyAssertion 时不要输出残缺对象，把候选 necessity
设为 uncertain、assertion_ids 设为空列表并在 unresolved_questions 说明缺口。
"""


async def load_bio_database_researcher_tools(
    config: ToolUniverseConfig | None = None,
) -> list[BaseTool]:
    """Load and verify the complete ToolUniverse allowlist."""

    available_tools = await load_tooluniverse_tools(config)
    tools_by_name = {tool.name: tool for tool in available_tools}
    missing = [
        name
        for name in BIO_DATABASE_RESEARCHER_TOOL_NAMES
        if name not in tools_by_name
    ]
    if missing:
        missing_names = ", ".join(missing)
        raise RuntimeError(
            f"ToolUniverse did not expose required tools: {missing_names}"
        )

    return [
        tools_by_name[name] for name in BIO_DATABASE_RESEARCHER_TOOL_NAMES
    ]


async def build_bio_database_researcher_subagent(
    *,
    config: ToolUniverseConfig | None = None,
    model: str | BaseChatModel | None = None,
    tools: Sequence[BaseTool] | None = None,
    middleware: Sequence[AgentMiddleware] | None = None,
) -> SubAgent:
    """Build the isolated Deep Agents database-research subagent."""

    subagent: SubAgent = {
        "name": BIO_DATABASE_RESEARCHER_NAME,
        "description": (
            "查询权威生物数据库，重建异源蛋白的原生催化系统、蛋白依赖和证据链。"
        ),
        "system_prompt": BIO_DATABASE_RESEARCHER_PROMPT,
        "tools": (
            list(tools)
            if tools is not None
            else await load_bio_database_researcher_tools(config)
        ),
        "response_format": BioDatabaseResearchResult,
    }
    if model is not None:
        subagent["model"] = model
    if middleware:
        subagent["middleware"] = list(middleware)
    return subagent
