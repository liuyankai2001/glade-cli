"""底盘感知的 KEGG 反向通路搜索工具。

本模块从目标 KEGG 化合物出发，沿“生成该化合物的反应”向上游反向展开，
直到所有主要前体都落入底盘 GEM 预先计算出的可达化合物集合。搜索过程中会：

1. 用 GEM 反应注释和通量上下界判断反应方向是否属于底盘内源能力；
2. 用 KEGG module 构造优先搜索的反应子图；
3. 按方向启发式、辅因子负担、电子风险和前体距离排列候选反应；
4. 用最大步数、最大异源酶数、状态支配和环路检测剪枝；
5. 将候选路线及逐步证据写入 CSV/JSON，交给后续 GEM 验证和酶系统选择。

注意：这里的搜索是启发式路线发现，不等同于热力学证明或最终通量可行性证明；
候选路线仍需经过下游 FBA/pFBA/FVA 验证。
"""

from __future__ import annotations

import heapq
import itertools
import json
import re
import time
from collections import OrderedDict
from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen
from langchain.tools import tool
from src.config.pathway_config import (
    DEFAULT_GAP_SEARCH_CONFIG,
    MAX_HETEROLOGOUS_ENZYMES,
    PATHWAY_RUNTIME_CONFIG,
)
from src.config.service_config import KEGG_HTTP_CONFIG, KEGG_REST_BASE_URL
from src.runtime.monitor import monitor
from src.tools.common.endogenous_directionality import (
    EndogenousDirectionIndex,
    direction_capable_reaction_count,
    is_directionally_endogenous,
    load_endogenous_direction_index,
)
from src.tools.common.session_paths import cache_dir as resolve_cache_dir
from src.tools.common.electron_transfer import (
    DEFAULT_ELECTRON_FALLBACK_MAX_NEW_ENZYMES,
    DEFAULT_ELECTRON_FALLBACK_MAX_REACTIONS_PER_COMPOUND,
    DEFAULT_ELECTRON_FALLBACK_MAX_ROUTES_PER_STATE,
    DEFAULT_ELECTRON_FALLBACK_MAX_TOTAL_STEPS,
    ELECTRON_AVOIDANCE_MODES,
    ElectronRequirement,
    HIGH_ELECTRON_RISK_SCORE,
    infer_electron_requirement,
    summarize_solution_electron_requirements,
)
from src.tools.common.session_paths import gem_model_file
from src.tools.common.session_paths import outputs_dir as resolve_outputs_dir
from src.tools.common.session_paths import producible_kegg_compounds_file
import cobra
import pandas as pd

from src.pathway_analyze.target_id import validate_target_compound_id


# ---------------------------------------------------------------------------
# 运行路径、网络请求和搜索默认值
# ---------------------------------------------------------------------------

CHASSIS_ANALYSIS_RUNNING_MARKER = ".analyze_chassis_metabolites.running"
REACHABLE_FILE_WAIT_SECONDS = PATHWAY_RUNTIME_CONFIG.reachable_file_wait_seconds
REACHABLE_FILE_POLL_SECONDS = PATHWAY_RUNTIME_CONFIG.reachable_file_poll_seconds
REACHABLE_FILE_STARTUP_GRACE_SECONDS = PATHWAY_RUNTIME_CONFIG.reachable_file_startup_grace_seconds
REACHABLE_MARKER_POLL_SECONDS = PATHWAY_RUNTIME_CONFIG.reachable_marker_poll_seconds

HTTP_TIMEOUT = KEGG_HTTP_CONFIG.timeout_seconds
HTTP_RETRIES = KEGG_HTTP_CONFIG.retries
REQUEST_SLEEP_SECONDS = KEGG_HTTP_CONFIG.sleep_seconds

DEFAULT_MAX_ROUTES_PER_STATE = DEFAULT_GAP_SEARCH_CONFIG.max_routes_per_state
DEFAULT_MAX_MODULE_CHAIN_DEPTH = DEFAULT_GAP_SEARCH_CONFIG.max_module_chain_depth
DEFAULT_MODULE_FILTER_MODE = DEFAULT_GAP_SEARCH_CONFIG.module_filter_mode
DEFAULT_ELECTRON_AVOIDANCE_MODE = DEFAULT_GAP_SEARCH_CONFIG.electron_avoidance_mode
SCREENING_RULE_VERSION = "directional_screening_v5_gem_bounds_cycle_safe"
REACTION_RESOLUTION_VERSION = "reaction_resolution.v1"
REACTION_RESOLUTION_MODES = frozenset({"strict", "audit"})
ENDOGENOUS_DIRECTION_MODE_GEM_BOUNDS = "gem_bounds"
CYCLE_POLICY = "no_reaction_reuse_no_expanded_compound_reentry"
ENDOGENOUS_DIRECTION_COMPARTMENTS = frozenset({"c"})  # 当前仅把胞质方向能力视为搜索锚点。

# 用于反应方向筛查和负担计算的 KEGG 通用化合物 ID。
OXYGEN_COMPOUND_ID = "C00007"
ATP_COMPOUND_ID = "C00002"
AMP_COMPOUND_ID = "C00020"
PPI_COMPOUND_ID = "C00013"
CO2_COMPOUND_ID = "C00011"
NADPH_COMPOUND_ID = "C00005"
SAM_COMPOUND_ID = "C00019"
SAH_COMPOUND_ID = "C00021"
COA_COMPOUND_ID = "C00010"
STEP_SUMMARY_COMMENT_PATTERN = re.compile(
    r"\b(?:multi|one|two|three|four|five|six|seven|eight|nine|ten|\d+)[ -]?steps?\b"
)
COMPONENT_STEP_COMMENT_PATTERN = re.compile(
    r"\b(?:the\s+)?(?:first|second|third|fourth|fifth|sixth|seventh|eighth|ninth|tenth|last)"
    r"\s+step\s+of\s+(?:the\s+)?(?:multi|one|two|three|four|five|six|seven|eight|nine|ten|\d+)"
    r"[ -]?step reaction\b"
    r"|\b(?:a\s+)?part\s+of\s+(?:the\s+)?(?:multi|one|two|three|four|five|six|seven|eight|nine|ten|\d+)"
    r"[ -]?step reaction\b"
)

# 这些通用小分子连接度极高。它们仍参与辅因子/电子风险计分，但不会作为
# “必须继续向上游解决的主要前体”，否则反向搜索树会被水、ATP、NADPH 等淹没。
IGNORED_COMMON_COMPOUNDS = {
    "C00001",  # H2O
    "C00002",  # ATP
    "C00003",  # NAD+
    "C00004",  # NADH
    "C00005",  # NADPH
    "C00006",  # NADP+
    "C00007",  # O2
    "C00008",  # ADP
    "C00009",  # Orthophosphate
    "C00010",  # CoA
    "C00011",  # CO2
    "C00013",  # PPi
    "C00014",  # NH3/NH4+
    "C00019",  # SAM
    "C00020",  # AMP
    "C00021",  # SAH
    "C00024",  # Acetyl-CoA
    "C00025",  # L-Glutamate
    "C00027",  # H2O2
    "C00028",  # Acceptor
    "C00030",  # Reduced acceptor
    "C00041",  # Alanine
    "C00044",  # GTP
    "C00063",  # CTP
    "C00064",  # L-Glutamine
    "C00075",  # UTP
    "C00080",  # H+
    "C00081",  # Sulfate
    "C00122",  # Fumarate
    "C00138",  # Reduced ferredoxin
    "C00139",  # Oxidized ferredoxin
    "C00229",  # ACP
    "C00268",  # Dihydrobiopterin
    "C00272",  # Tetrahydrobiopterin
    "C00342",  # Thioredoxin
    "C00343",  # Oxidized thioredoxin
    "C03024",  # Reduced NADPH---hemoprotein reductase
    "C03161",  # Oxidized NADPH---hemoprotein reductase
    "C01352",  # bicarbonate
    "C14818",  # reduced ferredoxin [iron-sulfur cluster]
}

KEGG_ID_PATTERN = re.compile(r"^[CDG]\d{5}$")
COMPOUND_TOKEN_PATTERN = re.compile(r"[CDG]\d{5}")
REACTION_ID_PATTERN = re.compile(r"R\d{5}")
MODULE_ID_PATTERN = re.compile(r"M\d{5}")


# ---------------------------------------------------------------------------
# 搜索过程中使用的不可变数据结构
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class ReactionRecord:
    """从 KEGG reaction 条目解析出的标准反应记录。"""

    reaction_id: str
    name: str
    comment: str
    equation: str
    annotation_text: str
    left_stoichiometry: Tuple[Tuple[str, float], ...]
    right_stoichiometry: Tuple[Tuple[str, float], ...]
    enzyme_ecs: Tuple[str, ...]
    ko_ids: Tuple[str, ...]
    pathway_ids: Tuple[str, ...]
    module_ids: Tuple[str, ...]
    rhea_ids: Tuple[str, ...] = tuple()


@dataclass(frozen=True)
class CompoundRecord:
    """KEGG compound 的名称及其通路、模块索引。"""

    name: str
    pathway_ids: Tuple[str, ...]
    module_ids: Tuple[str, ...]


@dataclass(frozen=True)
class ModuleRecord:
    """KEGG module 及其内部反应边、入口化合物和可生成化合物。"""

    name: str
    reaction_ids: Tuple[str, ...]
    reaction_edges: Tuple["ModuleReactionEdge", ...]
    start_compound_ids: Tuple[str, ...]
    produced_compound_ids: Tuple[str, ...]


@dataclass(frozen=True)
class ModuleReactionEdge:
    """KEGG module 中的一条有向反应边；一个边可能引用多个等价反应。"""

    reaction_ids: Tuple[str, ...]
    consumed_compound_ids: Tuple[str, ...]
    produced_compound_ids: Tuple[str, ...]


@dataclass(frozen=True)
class ReactionScreening:
    """一个反应方向的辅因子负担和启发式筛查结果。"""

    oxygen_required: bool
    nadph_burden: float
    sam_burden: float
    coa_burden: float
    thermo_direction: str
    rule_hits: Tuple[str, ...]


@dataclass(frozen=True)
class ReactionOption:
    """某反应以指定方向生成当前目标化合物时形成的一个搜索选项。

    ``precursor_compounds`` 是必须继续追溯的主要前体。水、辅因子等通用小分子
    在构造选项时被排除；筛查和电子需求各自收拢为一个不可变对象。
    """

    produced_compound: str
    reaction: ReactionRecord
    direction: str
    precursor_compounds: Tuple[str, ...]
    screening: ReactionScreening
    electron_requirement: ElectronRequirement


@dataclass(frozen=True)
class PlanStep:
    """已经加入候选路线的一步反应快照。

    静态反应信息由 ``option`` 复用；这里只保存路线相关的前体可达性和
    内源/异源状态。步骤编号由 ``Solution.steps`` 的顺序生成。
    """

    option: ReactionOption
    reachable_precursors: Tuple[str, ...]
    is_endogenous: bool
    source_reaction_ids: Tuple[str, ...] = tuple()
    resolution_action: str = "none"
    resolution_evidence: Tuple[str, ...] = tuple()


@dataclass(frozen=True)
class ReactionResolution:
    """KEGG 反应条目作为工程酶促步骤时的规范化状态。"""

    is_multistep: bool
    is_incomplete: bool
    component_ids: Tuple[str, ...]
    hard_blocker: bool
    reason: str


@dataclass(frozen=True)
class Solution:
    """一条所有主要前体均已连接到底盘可达集合的完整候选路线。"""

    steps: Tuple[PlanStep, ...]
    reaction_resolution_status: str = "resolved"
    normalization_events: Tuple[str, ...] = tuple()
    blocking_reaction_ids: Tuple[str, ...] = tuple()
    blocking_reasons: Tuple[str, ...] = tuple()

    @property
    def total_steps(self) -> int:
        """返回路线步骤总数。"""

        return len(self.steps)

    @property
    def heterologous_steps(self) -> int:
        """返回需要新增酶的步骤数。"""

        return sum(not step.is_endogenous for step in self.steps)

    @property
    def reaction_ready(self) -> bool:
        return not self.blocking_reaction_ids


@dataclass(frozen=True)
class SearchRoundResult:
    """单轮搜索产生的可推荐路线及被反应规范化门禁拒绝的路线。"""

    solutions: Tuple[Solution, ...]
    rejected_solutions: Tuple[Solution, ...]

    def __iter__(self):
        """兼容旧调用方把单轮结果直接当作 solution tuple 迭代。"""

        return iter(self.solutions)

    def __len__(self) -> int:
        return len(self.solutions)

    def __bool__(self) -> bool:
        return bool(self.solutions)


@dataclass(frozen=True, order=True)
class SearchQueueItem:
    """优先队列节点。

    dataclass 的字段顺序就是堆排序优先级：先比较异源反应数，再比较
    ``已走步数 + 未解决前体数``和未解决前体数。``sequence``仅用于
    稳定打破平局；总步数、具体状态和路径不参与对象比较。
    """

    heterologous_steps: int
    estimated_total_work: int
    unresolved_count: int
    sequence: int
    total_steps: int = field(compare=False)
    unresolved: frozenset[str] = field(compare=False)
    plan_steps: Tuple[PlanStep, ...] = field(compare=False)


@dataclass(frozen=True)
class SearchExecutionResult:
    """搜索结果及 module/电子风险兜底过程的可审计元数据。"""

    solutions: Tuple[Solution, ...]
    search_mode_used: str
    did_fallback: bool
    electron_avoidance_mode: str
    electron_avoidance_fallback: bool
    electron_fallback_parameters: Dict[str, int]
    target_compound_module_ids: Tuple[str, ...]
    matched_target_module_ids: Tuple[str, ...]
    direct_producing_reaction_ids: Tuple[str, ...]
    module_chain_context: "ModuleChainContext"
    rejected_solutions: Tuple[Solution, ...] = tuple()
    reaction_resolution_mode: str = "strict"


@dataclass(frozen=True)
class ModuleChainContext:
    """以目标 module 为起点向上游扩展得到的受限搜索子图。"""

    expanded_module_ids: Tuple[str, ...]
    boundary_compound_ids: Tuple[str, ...]
    allowed_reaction_ids: Tuple[str, ...]


# 路线签名由“本步产物、反应 ID、方向”组成，用于识别重复历史。
StateRouteSignature = Tuple[Tuple[str, str, str], ...]
# 状态记录保存“最大电子风险、异源步数、总步数、路线签名”，用于同状态
# Top-K 剪枝。电子风险必须进入状态成本：否则一条更短但依赖通用/未建模电子
# 载体的路线会把稍长、却可由 GEM 表达的显式低风险路线挤出搜索空间。
StateRouteRecord = Tuple[int, int, int, StateRouteSignature]


# ---------------------------------------------------------------------------
# 通用解析与 KEGG REST 客户端
# ---------------------------------------------------------------------------

def safe_mkdir(path: str | Path) -> Path:
    """确保目录存在并返回 Path，供缓存目录和输出目录统一使用。"""

    out = Path(path)
    out.mkdir(parents=True, exist_ok=True)
    return out


def normalize_annotation_value(value: Any) -> List[str]:
    """把 COBRA annotation 中可能出现的标量或集合统一转换为字符串列表。"""

    if value is None:
        return []
    if isinstance(value, str):
        return [value]
    if isinstance(value, (list, tuple, set)):
        return [str(item) for item in value if item is not None]
    return [str(value)]


def normalize_kegg_id(token: str) -> str:
    """去掉 ``cpd:``、``rn:`` 等命名空间前缀，保留标准 KEGG ID。"""

    token = token.strip()
    if ":" in token:
        token = token.split(":", 1)[1]
    return token


def kegg_compound_namespace(compound_id: str) -> str:
    """根据 C/D/G 前缀选择 KEGG compound/drug/glycan 命名空间。"""

    compound_id = normalize_kegg_id(compound_id)
    if compound_id.startswith("G"):
        return "gl"
    if compound_id.startswith("D"):
        return "dr"
    return "cpd"


def dedupe_keep_order(values: Iterable[str]) -> Tuple[str, ...]:
    """稳定去重：保留第一次出现的顺序，保证搜索和输出可复现。"""

    seen = set()
    ordered: List[str] = []
    for value in values:
        if value in seen:
            continue
        seen.add(value)
        ordered.append(value)
    return tuple(ordered)


def parse_kegg_flatfile(text: str) -> Dict[str, List[str]]:
    """解析 KEGG 固定 12 字符字段宽度的文本格式。

    同一字段的续行会追加到列表中，直到遇到条目终止符 ``///``。
    """

    fields: Dict[str, List[str]] = {}
    current_key: str | None = None

    for raw_line in text.splitlines():
        if raw_line.strip() == "///":
            break

        key = raw_line[:12].strip()
        value = raw_line[12:].rstrip()

        if key:
            current_key = key
            fields.setdefault(key, []).append(value.strip())
        elif current_key is not None:
            fields[current_key].append(value.strip())

    return fields


def first_or_empty(values: Sequence[str]) -> str:
    """返回序列首项；空序列统一返回空字符串。"""

    return values[0] if values else ""


class KeggRestClient:
    """带内存缓存、磁盘缓存、批量预取和重试的轻量 KEGG REST 客户端。

    默认使用类级共享缓存，使一次任务中创建的多个客户端可以复用条目；测试可
    通过 ``use_shared_cache=False`` 隔离状态。磁盘缓存用于减少重复网络请求。
    """

    _GLOBAL_TEXT_CACHE: Dict[str, str] = {}
    _GLOBAL_REACTION_CACHE: Dict[str, ReactionRecord] = {}
    _GLOBAL_COMPOUND_RECORD_CACHE: Dict[str, CompoundRecord] = {}
    _GLOBAL_MODULE_CACHE: Dict[str, ModuleRecord] = {}
    _GLOBAL_COMPOUND_REACTIONS_CACHE: Dict[str, Tuple[str, ...]] = {}

    def __init__(
        self,
        cache_dir: str | Path | None = None,
        timeout: int = HTTP_TIMEOUT,
        request_sleep_seconds: float = REQUEST_SLEEP_SECONDS,
        use_shared_cache: bool = True,
    ) -> None:
        """配置缓存目录、请求超时和限速，并选择共享或实例级缓存。"""

        self.cache_dir = safe_mkdir(cache_dir).resolve() if cache_dir is not None else None
        self.timeout = timeout
        self.request_sleep_seconds = max(0.0, float(request_sleep_seconds))

        if use_shared_cache:
            self._text_cache = self._GLOBAL_TEXT_CACHE
            self._reaction_cache = self._GLOBAL_REACTION_CACHE
            self._compound_record_cache = self._GLOBAL_COMPOUND_RECORD_CACHE
            self._module_cache = self._GLOBAL_MODULE_CACHE
            self._compound_reactions_cache = self._GLOBAL_COMPOUND_REACTIONS_CACHE
        else:
            self._text_cache = {}
            self._reaction_cache = {}
            self._compound_record_cache = {}
            self._module_cache = {}
            self._compound_reactions_cache = {}

    def _cache_path(self, namespace: str, key: str, suffix: str) -> Path | None:
        """构造单个缓存文件路径；未启用磁盘缓存时返回 ``None``。"""

        if self.cache_dir is None:
            return None
        return safe_mkdir(self.cache_dir / namespace) / f"{key}{suffix}"

    def _fetch_text(self, url: str, cache_key: str) -> str:
        """按“内存 → 磁盘 → KEGG 网络”的顺序读取文本，并在网络失败时重试。"""

        if cache_key in self._text_cache:
            return self._text_cache[cache_key]

        namespace, raw_key = cache_key.split(":", 1)
        cache_path = self._cache_path(namespace, raw_key, ".txt")
        if cache_path is not None and cache_path.exists():
            text = cache_path.read_text(encoding="utf-8")
            self._text_cache[cache_key] = text
            return text

        last_error: Exception | None = None
        for attempt in range(1, HTTP_RETRIES + 1):
            try:
                request = Request(url, headers={"User-Agent": "Mozilla/5.0"})
                with urlopen(request, timeout=self.timeout) as response:
                    text = response.read().decode("utf-8")
                self._text_cache[cache_key] = text
                if cache_path is not None:
                    cache_path.write_text(text, encoding="utf-8")
                if self.request_sleep_seconds > 0:
                    time.sleep(self.request_sleep_seconds)
                return text
            except (HTTPError, URLError, TimeoutError, ConnectionError) as exc:
                last_error = exc
                time.sleep(attempt)

        raise RuntimeError(f"Failed to fetch KEGG URL after retries: {url}") from last_error

    def get_compound_name(self, compound_id: str) -> str:
        """返回化合物主名称；化合物记录缓存会避免重复请求。"""

        return self.get_compound_record(compound_id).name

    def get_compound_record(self, compound_id: str) -> CompoundRecord:
        """下载并解析化合物条目，提取 pathway、module 和关联反应。"""

        compound_id = normalize_kegg_id(compound_id)
        if compound_id in self._compound_record_cache:
            return self._compound_record_cache[compound_id]

        namespace = kegg_compound_namespace(compound_id)
        text = self._fetch_text(
            f"{KEGG_REST_BASE_URL}/get/{namespace}:{compound_id}",
            f"compound:{compound_id}",
        )
        fields = parse_kegg_flatfile(text)

        name = first_or_empty(fields.get("NAME", []))
        if ";" in name:
            name = name.split(";", 1)[0]

        record = CompoundRecord(
            name=name,
            pathway_ids=parse_pathway_field(fields.get("PATHWAY", [])),
            module_ids=parse_module_field(fields.get("MODULE", [])),
        )
        self._compound_record_cache[compound_id] = record
        return record

    def get_reaction_ids_for_compound(self, compound_id: str) -> Tuple[str, ...]:
        """调用 KEGG link 接口获取所有包含该化合物的反应 ID。"""

        compound_id = normalize_kegg_id(compound_id)
        if compound_id in self._compound_reactions_cache:
            return self._compound_reactions_cache[compound_id]

        namespace = kegg_compound_namespace(compound_id)
        text = self._fetch_text(
            f"{KEGG_REST_BASE_URL}/link/reaction/{namespace}:{compound_id}",
            f"compound_to_reaction:{compound_id}",
        )

        reaction_ids: List[str] = []
        for line in text.splitlines():
            if "\t" not in line:
                continue
            _, raw_reaction_id = line.split("\t", 1)
            reaction_ids.append(normalize_kegg_id(raw_reaction_id))

        ordered = dedupe_keep_order(reaction_ids)
        self._compound_reactions_cache[compound_id] = ordered
        return ordered

    def get_reaction(self, reaction_id: str) -> ReactionRecord:
        """获取并解析一个 KEGG reaction；已缓存条目不会再次请求网络。"""

        reaction_id = normalize_kegg_id(reaction_id)
        if reaction_id in self._reaction_cache:
            return self._reaction_cache[reaction_id]

        text = self._fetch_text(
            f"{KEGG_REST_BASE_URL}/get/rn:{reaction_id}",
            f"reaction:{reaction_id}",
        )
        record = self._parse_reaction_text(reaction_id, text)
        self._reaction_cache[reaction_id] = record
        return record

    def try_get_reaction(self, reaction_id: str) -> ReactionRecord | None:
        """容错获取反应；网络或解析失败时返回 ``None`` 供调用方跳过。"""

        try:
            return self.get_reaction(reaction_id)
        except RuntimeError:
            return None

    def prefetch_reactions(self, reaction_ids: Sequence[str], batch_size: int = 10) -> None:
        """批量预取反应，批量请求失败时退回逐条请求，避免整批数据丢失。"""

        missing = [
            normalize_kegg_id(reaction_id)
            for reaction_id in reaction_ids
            if normalize_kegg_id(reaction_id) not in self._reaction_cache
        ]
        if not missing:
            return

        for start in range(0, len(missing), batch_size):
            batch = missing[start : start + batch_size]
            url_ids = "+".join(f"rn:{reaction_id}" for reaction_id in batch)
            try:
                batch_text = self._fetch_text(
                    f"{KEGG_REST_BASE_URL}/get/{url_ids}",
                    f"reaction_batch:{'_'.join(batch)}",
                )
            except RuntimeError:
                for reaction_id in batch:
                    self.try_get_reaction(reaction_id)
                continue
            for entry_text in split_kegg_entries(batch_text):
                entry_fields = parse_kegg_flatfile(entry_text)
                entry_id = first_or_empty(entry_fields.get("ENTRY", [])).split()[0]
                if not entry_id:
                    continue
                normalized_entry_id = normalize_kegg_id(entry_id)
                self._text_cache[f"reaction:{normalized_entry_id}"] = entry_text
                reaction_cache_path = self._cache_path("reaction", normalized_entry_id, ".txt")
                if reaction_cache_path is not None and not reaction_cache_path.exists():
                    reaction_cache_path.write_text(entry_text, encoding="utf-8")
                self._reaction_cache[normalized_entry_id] = self._parse_reaction_text(
                    normalized_entry_id,
                    entry_text,
                )

            for reaction_id in batch:
                if reaction_id in self._reaction_cache:
                    continue
                self.try_get_reaction(reaction_id)

    def get_module(self, module_id: str) -> ModuleRecord:
        """解析 KEGG module，并从 REACTION 字段恢复其有向反应图。"""

        module_id = normalize_kegg_id(module_id)
        if module_id in self._module_cache:
            return self._module_cache[module_id]

        text = self._fetch_text(
            f"{KEGG_REST_BASE_URL}/get/md:{module_id}",
            f"module:{module_id}",
        )
        fields = parse_kegg_flatfile(text)
        reaction_edges = parse_module_reaction_edges(fields.get("REACTION", []))
        module_reaction_ids = dedupe_keep_order(
            reaction_id
            for edge in reaction_edges
            for reaction_id in edge.reaction_ids
        )
        produced_compound_ids = dedupe_keep_order(
            compound_id
            for edge in reaction_edges
            for compound_id in edge.produced_compound_ids
        )
        consumed_compound_ids = dedupe_keep_order(
            compound_id
            for edge in reaction_edges
            for compound_id in edge.consumed_compound_ids
        )
        produced_compound_set = set(produced_compound_ids)

        record = ModuleRecord(
            name=first_or_empty(fields.get("NAME", [])),
            reaction_ids=module_reaction_ids,
            reaction_edges=reaction_edges,
            start_compound_ids=dedupe_keep_order(
                compound_id
                for compound_id in consumed_compound_ids
                if compound_id not in produced_compound_set
            ),
            produced_compound_ids=produced_compound_ids,
        )
        self._module_cache[module_id] = record
        return record

    def prefetch_modules(self, module_ids: Sequence[str]) -> None:
        """逐个预取 module；底层缓存会自动消除重复请求。"""

        for module_id in dedupe_keep_order(normalize_kegg_id(module_id) for module_id in module_ids):
            self.get_module(module_id)

    def _parse_reaction_text(self, reaction_id: str, text: str) -> ReactionRecord:
        """把一个 KEGG reaction 平面文本转换为标准 ``ReactionRecord``。"""

        fields = parse_kegg_flatfile(text)

        equation = " ".join(fields.get("EQUATION", []))
        left_side, right_side = split_equation(equation)
        left_stoichiometry = parse_equation_side(left_side)
        right_stoichiometry = parse_equation_side(right_side)
        enzyme_ecs = parse_enzyme_field(fields.get("ENZYME", []))
        ko_ids = parse_orthology_field(fields.get("ORTHOLOGY", []))
        pathway_ids = parse_pathway_field(fields.get("PATHWAY", []))
        module_ids = parse_module_field(fields.get("MODULE", []))
        rhea_ids = parse_rhea_field(fields.get("DBLINKS", []))
        annotation_text = " ".join(
            " ".join(fields.get(field_name, []))
            for field_name in ("DEFINITION", "PATHWAY", "MODULE", "ORTHOLOGY")
        )

        record = ReactionRecord(
            reaction_id=reaction_id,
            name=first_or_empty(fields.get("NAME", [])),
            comment=" ".join(fields.get("COMMENT", [])),
            equation=equation,
            annotation_text=annotation_text,
            left_stoichiometry=left_stoichiometry,
            right_stoichiometry=right_stoichiometry,
            enzyme_ecs=enzyme_ecs,
            ko_ids=ko_ids,
            pathway_ids=pathway_ids,
            module_ids=module_ids,
            rhea_ids=rhea_ids,
        )
        return record


# ---------------------------------------------------------------------------
# 反应方程式解析、方向启发式与辅因子负担
# ---------------------------------------------------------------------------

def split_equation(equation: str) -> Tuple[str, str]:
    """按 KEGG 箭头拆分方程左右两侧，不在这里改变反应方向。"""

    for arrow in ("<=>", "=>", "<="):
        if arrow in equation:
            left, right = equation.split(arrow, 1)
            return left.strip(), right.strip()
    return equation.strip(), ""


def parse_equation_side(side: str) -> Tuple[Tuple[str, float], ...]:
    """解析方程一侧的 KEGG 化合物及计量系数。

    KEGG 中无法解析或缺失的系数按 1 处理；相同化合物重复出现时合并系数。
    """

    amounts: "OrderedDict[str, float]" = OrderedDict()
    for raw_term in side.split("+"):
        term = raw_term.strip()
        if not term:
            continue
        match = COMPOUND_TOKEN_PATTERN.search(term)
        if not match:
            continue
        compound_id = match.group(0)
        coefficient = 1.0
        prefix = term[: match.start()].strip()
        if prefix:
            prefix_token = prefix.split()[-1]
            try:
                coefficient = float(prefix_token)
            except ValueError:
                coefficient = 1.0
        amounts[compound_id] = amounts.get(compound_id, 0.0) + coefficient
    return tuple((compound_id, value) for compound_id, value in amounts.items())


def stoichiometry_amount(stoichiometry: Sequence[Tuple[str, float]], compound_id: str) -> float:
    """返回指定化合物在一侧计量关系中的系数，缺失时为 0。"""

    for current_compound_id, amount in stoichiometry:
        if current_compound_id == compound_id:
            return amount
    return 0.0


def stoichiometry_compound_ids(
    stoichiometry: Sequence[Tuple[str, float]],
) -> Tuple[str, ...]:
    """从计量关系中提取化合物 ID。"""

    return tuple(compound_id for compound_id, _ in stoichiometry)


def lower_reaction_text(reaction: ReactionRecord) -> str:
    """拼接反应名称、定义和注释，形成小写关键词检索文本。"""

    return " ".join(
        part
        for part in (
            reaction.name,
            reaction.comment,
            reaction.annotation_text,
        )
        if part
    ).lower()


def left_to_right_signature(reaction: ReactionRecord) -> Dict[str, bool]:
    """提取原始左→右方向的反应类型特征，用于方向合理性启发式。

    这些特征只基于方程和文字标签，并不是标准反应自由能计算。
    """

    text = lower_reaction_text(reaction)
    left_consumed = {compound_id for compound_id, amount in reaction.left_stoichiometry if amount > 0}
    right_produced = {compound_id for compound_id, amount in reaction.right_stoichiometry if amount > 0}
    consumes_atp = ATP_COMPOUND_ID in left_consumed
    consumes_coa = COA_COMPOUND_ID in left_consumed
    consumes_sam = SAM_COMPOUND_ID in left_consumed
    consumes_oxygen = OXYGEN_COMPOUND_ID in left_consumed
    produces_amp_or_ppi = AMP_COMPOUND_ID in right_produced or PPI_COMPOUND_ID in right_produced
    produces_sah = SAH_COMPOUND_ID in right_produced
    produces_co2 = CO2_COMPOUND_ID in right_produced

    oxygen_signature = any(
        token in text
        for token in (
            "oxygen oxidoreductase",
            "oxygenase",
            "hydroxylase",
            "monooxygenase",
            "dioxygenase",
        )
    )
    methyltransferase_signature = "methyltransferase" in text or (consumes_sam and produces_sah)
    amp_forming_signature = "amp-forming" in text or (consumes_atp and produces_amp_or_ppi)
    coa_activation_signature = (
        "coa ligase" in text
        or "coenzyme a ligase" in text
        or ("coa-acylating" in text)
        or (consumes_coa and consumes_atp and produces_amp_or_ppi)
    )
    decarboxylating_only_signature = "decarboxylating" in text and "carboxylating" not in text
    atp_ppi_coupling_signature = consumes_atp or produces_amp_or_ppi

    return {
        "oxygen_signature": oxygen_signature,
        "methyltransferase_signature": methyltransferase_signature,
        "amp_forming_signature": amp_forming_signature,
        "coa_activation_signature": coa_activation_signature,
        "decarboxylating_only_signature": decarboxylating_only_signature,
        "atp_ppi_coupling_signature": atp_ppi_coupling_signature,
        "consumes_oxygen": consumes_oxygen,
        "consumes_sam": consumes_sam,
        "consumes_coa": consumes_coa,
        "produces_co2": produces_co2,
    }


def directional_stoichiometry(
    reaction: ReactionRecord,
    direction: str,
) -> Tuple[Tuple[Tuple[str, float], ...], Tuple[Tuple[str, float], ...]]:
    """按照当前搜索方向返回“消耗侧、生成侧”的计量关系。"""

    if direction == "left_to_right":
        return reaction.left_stoichiometry, reaction.right_stoichiometry
    if direction == "right_to_left":
        return reaction.right_stoichiometry, reaction.left_stoichiometry
    raise ValueError(f"Unsupported reaction direction: {direction}")


def burden_value(stoichiometry: Sequence[Tuple[str, float]], compound_id: str) -> float:
    """读取当前消耗侧某辅因子的计量负担。"""

    return stoichiometry_amount(stoichiometry, compound_id)


def screening_hits_to_text(rule_hits: Sequence[str]) -> str:
    """把筛查规则命中项序列化为适合 CSV 展示的文本。"""

    return ";".join(rule_hits)


def screen_reaction_direction(
    *,
    reaction: ReactionRecord,
    direction: str,
    consumed_stoichiometry: Tuple[Tuple[str, float], ...],
) -> Tuple[bool, float, float, float, str, Tuple[str, ...]]:
    """评估一个具体反应方向的耗氧、辅因子负担和方向启发式。

    返回值依次为：是否耗氧、NADPH 消耗、SAM 消耗、CoA 消耗、
    方向启发式标签以及命中的可审计规则。``thermo_direction`` 只是根据
    氧化、脱羧、ATP/PPi 偶联等模式得到的 favored/neutral/disfavored，
    不能替代 ΔG 或热力学可行性计算。
    """

    rule_hits: List[str] = []
    signature = left_to_right_signature(reaction)
    oxygen_required = burden_value(consumed_stoichiometry, OXYGEN_COMPOUND_ID) > 0
    if oxygen_required:
        rule_hits.append("oxygen_required")

    # 只统计当前方向实际位于消耗侧的计量系数，不把生成的辅因子记作负担。
    nadph_burden = burden_value(consumed_stoichiometry, NADPH_COMPOUND_ID)
    sam_burden = burden_value(consumed_stoichiometry, SAM_COMPOUND_ID)
    coa_burden = burden_value(consumed_stoichiometry, COA_COMPOUND_ID)

    # 对典型强方向反应进行反向使用时留下风险证据，便于结果审计。
    if direction == "right_to_left":
        if signature["oxygen_signature"]:
            rule_hits.append("reverse_oxygenase_signature")
        if signature["amp_forming_signature"] or signature["coa_activation_signature"]:
            rule_hits.append("reverse_amp_forming_or_coa_activation_signature")
        if signature["methyltransferase_signature"]:
            rule_hits.append("reverse_methyltransferase_signature")
        if signature["decarboxylating_only_signature"]:
            rule_hits.append("reverse_decarboxylating_only_signature")

    thermo_hits: List[str] = []
    if direction == "left_to_right":
        if signature["oxygen_signature"] and (signature["consumes_oxygen"] or oxygen_required):
            thermo_hits.append("thermo_favored_forward_oxygenation")
        if signature["decarboxylating_only_signature"] and signature["produces_co2"]:
            thermo_hits.append("thermo_favored_forward_decarboxylation")
        if signature["amp_forming_signature"] or signature["atp_ppi_coupling_signature"]:
            thermo_hits.append("thermo_favored_forward_atp_ppi_coupling")
        if signature["methyltransferase_signature"] and signature["consumes_sam"]:
            thermo_hits.append("thermo_favored_forward_sam_transfer")
        if signature["coa_activation_signature"] and signature["consumes_coa"]:
            thermo_hits.append("thermo_favored_forward_coa_activation")
    else:
        if signature["atp_ppi_coupling_signature"]:
            thermo_hits.append("thermo_disfavored_reverse_atp_ppi_coupling")
        if signature["methyltransferase_signature"]:
            thermo_hits.append("thermo_disfavored_reverse_sam_transfer")
        if signature["coa_activation_signature"]:
            thermo_hits.append("thermo_disfavored_reverse_coa_activation")
        if signature["decarboxylating_only_signature"]:
            thermo_hits.append("thermo_disfavored_reverse_decarboxylation")

    # 同一方向可能命中多条规则；存在 favored 证据时优先标为 favored。
    thermo_direction = "neutral"
    if any(hit.startswith("thermo_favored") for hit in thermo_hits):
        thermo_direction = "favored"
    elif any(hit.startswith("thermo_disfavored") for hit in thermo_hits):
        thermo_direction = "disfavored"

    rule_hits.extend(thermo_hits)

    return (
        oxygen_required,
        nadph_burden,
        sam_burden,
        coa_burden,
        thermo_direction,
        tuple(rule_hits),
    )


# ---------------------------------------------------------------------------
# KEGG 字段提取与 module 反应图解析
# ---------------------------------------------------------------------------

def split_kegg_entries(text: str) -> List[str]:
    """把 KEGG 批量 ``get`` 响应按 ``///`` 拆成独立条目。"""

    entries: List[str] = []
    buffer: List[str] = []
    for line in text.splitlines():
        buffer.append(line)
        if line.strip() == "///":
            entries.append("\n".join(buffer) + "\n")
            buffer = []
    if buffer:
        entries.append("\n".join(buffer) + "\n")
    return entries


def parse_enzyme_field(lines: Sequence[str]) -> Tuple[str, ...]:
    """从 KEGG ENZYME 字段提取并稳定去重 EC 编号。"""

    values: List[str] = []
    for line in lines:
        values.extend(token for token in line.split() if re.match(r"^\d+\.\d+\..+", token))
    return dedupe_keep_order(values)


def parse_orthology_field(lines: Sequence[str]) -> Tuple[str, ...]:
    """从 ORTHOLOGY 字段提取 KO 编号。"""

    values: List[str] = []
    for line in lines:
        match = re.match(r"^(K\d{5})\b", line)
        if match:
            values.append(match.group(1))
    return dedupe_keep_order(values)


def parse_pathway_field(lines: Sequence[str]) -> Tuple[str, ...]:
    """从 PATHWAY 字段提取并规范化 KEGG pathway 编号。"""

    values: List[str] = []
    for line in lines:
        match = re.match(r"^(?:rn|map)(\d{5})\b", line)
        if match:
            values.append(f"map{match.group(1)}")
    return dedupe_keep_order(values)


def parse_module_field(lines: Sequence[str]) -> Tuple[str, ...]:
    """从 MODULE 字段提取 KEGG module 编号。"""

    values: List[str] = []
    for line in lines:
        values.extend(MODULE_ID_PATTERN.findall(line))
    return dedupe_keep_order(values)


def parse_rhea_field(lines: Sequence[str]) -> Tuple[str, ...]:
    """Extract normalized Rhea identifiers from KEGG ``DBLINKS`` lines."""

    values: List[str] = []
    for line in lines:
        match = re.match(r"^RHEA:\s*(.*)$", str(line).strip(), flags=re.IGNORECASE)
        if not match:
            continue
        values.extend(
            f"RHEA:{identifier}"
            for identifier in re.findall(r"\d+", match.group(1))
        )
    return dedupe_keep_order(values)


def extract_reaction_ids_from_text(text: str) -> Tuple[str, ...]:
    """从任意 KEGG 文本片段中提取反应编号。"""

    return dedupe_keep_order(REACTION_ID_PATTERN.findall(text))


def parse_module_reaction_edges(lines: Sequence[str]) -> Tuple[ModuleReactionEdge, ...]:
    """把 module 的 REACTION 行解析为有向的“消耗物→生成物”边。"""

    edges: List[ModuleReactionEdge] = []
    for line in lines:
        reaction_ids = extract_reaction_ids_from_text(line)
        compound_match = COMPOUND_TOKEN_PATTERN.search(line)
        if not reaction_ids or not compound_match:
            continue

        equation_text = line[compound_match.start() :].strip()
        if "<=>" in equation_text:
            left_side, right_side = equation_text.split("<=>", 1)
        elif "->" in equation_text:
            left_side, right_side = equation_text.split("->", 1)
        elif "=>" in equation_text:
            left_side, right_side = equation_text.split("=>", 1)
        else:
            continue

        consumed = dedupe_keep_order(compound_id for compound_id, _ in parse_equation_side(left_side))
        produced = dedupe_keep_order(compound_id for compound_id, _ in parse_equation_side(right_side))
        if not consumed and not produced:
            continue
        edges.append(
            ModuleReactionEdge(
                reaction_ids=reaction_ids,
                consumed_compound_ids=consumed,
                produced_compound_ids=produced,
            )
        )
    return tuple(edges)


def reaction_option_allowed_by_modules(
    option: ReactionOption,
    allowed_reaction_ids: set[str] | None,
) -> bool:
    """判断候选反应是否属于当前 module 链允许的反应集合。

    KEGG 的多步汇总反应本身可能不在 module 中；如果其注释引用的所有组件
    反应都在允许集合中，也允许该汇总选项进入后续拆解。
    """

    if allowed_reaction_ids is None:
        return True
    if not allowed_reaction_ids:
        return False
    if option.reaction.reaction_id in allowed_reaction_ids:
        return True

    referenced_reaction_ids = extract_reaction_ids_from_text(option.reaction.comment)
    return bool(referenced_reaction_ids) and all(
        reaction_id in allowed_reaction_ids
        for reaction_id in referenced_reaction_ids
    )


# ---------------------------------------------------------------------------
# 多步汇总反应拆解与单步 ReactionOption 构造
# ---------------------------------------------------------------------------

def multistep_component_reaction_ids(reaction: ReactionRecord) -> Tuple[str, ...]:
    """从 KEGG COMMENT 中提取多步汇总反应引用的真实组件反应 ID。"""

    comment = reaction.comment.lower()
    if not STEP_SUMMARY_COMMENT_PATTERN.search(comment):
        return tuple()
    component_ids = tuple(
        reaction_id
        for reaction_id in extract_reaction_ids_from_text(reaction.comment)
        if reaction_id != reaction.reaction_id
    )
    return component_ids


def is_component_step_comment(reaction: ReactionRecord) -> bool:
    """判断条目注释是否表明它已经是多步反应中的一个组件步骤。"""

    return bool(COMPONENT_STEP_COMMENT_PATTERN.search(reaction.comment.lower()))


def should_expand_multistep_summary(reaction: ReactionRecord) -> bool:
    """仅拆解汇总条目，避免把已经是组件步骤的条目再次递归拆解。"""

    if is_component_step_comment(reaction):
        return False
    return bool(multistep_component_reaction_ids(reaction))


def classify_reaction_resolution(reaction: ReactionRecord) -> ReactionResolution:
    """判断 KEGG 条目能否直接作为一个可实现的酶促步骤。

    ``incomplete reaction`` 且没有 EC/KO 表示无法锁定独立工程步骤；若条目
    已有明确 EC/KO，则保留为 enzyme-anchored review，交由严格反应证据门禁
    继续验证。仅标记 ``multistep`` 但具有 EC/KO 的条目也可能表示同一酶的
    多轮催化，不作硬阻断。
    """

    comment = reaction.comment.lower()
    is_multistep = bool(STEP_SUMMARY_COMMENT_PATTERN.search(comment))
    is_incomplete = "incomplete reaction" in comment
    component_ids = multistep_component_reaction_ids(reaction)
    unresolved_multistep = is_multistep and not component_ids
    lacks_enzyme_anchor = not reaction.enzyme_ecs and not reaction.ko_ids
    hard_blocker = (
        is_incomplete and lacks_enzyme_anchor
    ) or (unresolved_multistep and lacks_enzyme_anchor)
    if is_incomplete and lacks_enzyme_anchor:
        reason = "incomplete_reaction"
    elif is_incomplete:
        reason = "incomplete_enzyme_anchored_review"
    elif unresolved_multistep and lacks_enzyme_anchor:
        reason = "unresolved_multistep_without_enzyme_anchor"
    elif unresolved_multistep:
        reason = "multistep_enzyme_anchored_review"
    elif component_ids:
        reason = "explicit_multistep"
    else:
        reason = "atomic"
    return ReactionResolution(
        is_multistep=is_multistep,
        is_incomplete=is_incomplete,
        component_ids=component_ids,
        hard_blocker=hard_blocker,
        reason=reason,
    )


def make_reaction_option(
    compound_id: str,
    reaction: ReactionRecord,
    direction: str,
    ignored_common_compounds: set[str],
) -> ReactionOption:
    """把一个 KEGG 反应的确定方向转换为可排序、可扩展的搜索选项。"""

    # 反向搜索只关心“为了生成 compound_id，需要消耗哪些前体”。
    consumed_stoichiometry, produced_stoichiometry = directional_stoichiometry(reaction, direction)
    precursors = dedupe_keep_order(
        precursor_id
        for precursor_id, _ in consumed_stoichiometry
        if precursor_id != compound_id
    )
    # 通用辅因子不进入未解决状态，但会在下方继续计算负担和电子风险。
    tracked = dedupe_keep_order(x for x in precursors if x not in ignored_common_compounds)
    (
        oxygen_required,
        nadph_burden,
        sam_burden,
        coa_burden,
        thermo_direction,
        screening_rule_hits,
    ) = screen_reaction_direction(
        reaction=reaction,
        direction=direction,
        consumed_stoichiometry=consumed_stoichiometry,
    )
    # 电子模块只负责识别和标注风险，具体电子伙伴由下游酶系统选择解决。
    electron_requirement = infer_electron_requirement(
        consumed_stoichiometry=consumed_stoichiometry,
        produced_stoichiometry=produced_stoichiometry,
        reaction_name=reaction.name,
        reaction_comment=reaction.comment,
        equation=reaction.equation,
        enzyme_ecs=reaction.enzyme_ecs,
    )
    return ReactionOption(
        produced_compound=compound_id,
        reaction=reaction,
        direction=direction,
        precursor_compounds=tracked,
        screening=ReactionScreening(
            oxygen_required=oxygen_required,
            nadph_burden=nadph_burden,
            sam_burden=sam_burden,
            coa_burden=coa_burden,
            thermo_direction=thermo_direction,
            rule_hits=screening_rule_hits,
        ),
        electron_requirement=electron_requirement,
    )


def reaction_compounds(reaction: ReactionRecord) -> Tuple[str, ...]:
    """合并反应两侧出现的化合物并稳定去重。"""

    return dedupe_keep_order(
        compound_id
        for stoichiometry in (reaction.left_stoichiometry, reaction.right_stoichiometry)
        for compound_id in stoichiometry_compound_ids(stoichiometry)
    )


def infer_direction_for_product(reaction: ReactionRecord, product_compound: str) -> Optional[str]:
    """推断生成指定产物的唯一方向；两侧都含该产物时返回 ``None``。"""

    in_left = stoichiometry_amount(reaction.left_stoichiometry, product_compound) > 0
    in_right = stoichiometry_amount(reaction.right_stoichiometry, product_compound) > 0
    if in_left == in_right:
        return None
    return "right_to_left" if in_left else "left_to_right"


def add_screening_hits(option: ReactionOption, extra_hits: Sequence[str]) -> ReactionOption:
    """复制搜索选项并附加新的可审计筛查证据。"""

    return replace(
        option,
        screening=replace(
            option.screening,
            rule_hits=dedupe_keep_order((*option.screening.rule_hits, *extra_hits)),
        ),
    )


def fallback_multistep_option(option: ReactionOption) -> Tuple[ReactionOption, ...]:
    """为无法拆解的汇总反应添加标记并作为单步候选保留。"""

    return (
        add_screening_hits(
            option,
            (f"multistep_decomposition_fallback:{option.reaction.reaction_id}",),
        ),
    )


def candidate_multistep_orders(component_ids: Tuple[str, ...]) -> Tuple[Tuple[str, ...], ...]:
    """生成组件正序和逆序候选，避免假设注释顺序就是合成方向。"""

    orders: List[Tuple[str, ...]] = []
    for order in (tuple(reversed(component_ids)), component_ids):
        if order not in orders:
            orders.append(order)
    return tuple(orders)


def try_decompose_multistep_order(
    option: ReactionOption,
    component_ids: Tuple[str, ...],
    client: KeggRestClient,
    ignored_common_compounds: set[str],
) -> Tuple[ReactionOption, ...] | None:
    """尝试按给定顺序连接组件反应；任一中间体无法唯一衔接即判定失败。"""

    current_compound = option.produced_compound
    component_options: List[ReactionOption] = []

    for index, reaction_id in enumerate(component_ids):
        reaction = client.try_get_reaction(reaction_id)
        if reaction is None:
            return None

        direction = infer_direction_for_product(reaction, current_compound)
        if direction is None:
            return None

        component_option = make_reaction_option(
            compound_id=current_compound,
            reaction=reaction,
            direction=direction,
            ignored_common_compounds=ignored_common_compounds,
        )
        component_option = add_screening_hits(
            component_option,
            (
                f"decomposed_from:{option.reaction.reaction_id}",
                f"decomposition_index:{index + 1}/{len(component_ids)}",
            ),
        )
        component_options.append(component_option)

        if index >= len(component_ids) - 1:
            continue

        next_reaction = client.try_get_reaction(component_ids[index + 1])
        if next_reaction is None:
            return None
        next_compounds = set(reaction_compounds(next_reaction))
        next_intermediates = [
            compound_id
            for compound_id in component_option.precursor_compounds
            if compound_id in next_compounds
        ]
        if len(next_intermediates) != 1:
            return None
        current_compound = next_intermediates[0]

    return tuple(component_options)


def decompose_multistep_option(
    option: ReactionOption,
    client: KeggRestClient,
    ignored_common_compounds: set[str],
) -> Tuple[ReactionOption, ...]:
    """将 KEGG 多步汇总反应替换为逐步反应；无法可靠拆解时保留原条目。"""

    resolution = classify_reaction_resolution(option.reaction)
    component_ids = resolution.component_ids
    if resolution.is_multistep and not component_ids:
        return (
            add_screening_hits(
                option,
                (
                    f"reaction_resolution:{resolution.reason}",
                    f"multistep_decomposition_fallback:{option.reaction.reaction_id}",
                ),
            ),
        )
    if not should_expand_multistep_summary(option.reaction) or not component_ids:
        if resolution.is_incomplete:
            return (
                add_screening_hits(
                    option,
                    (f"reaction_resolution:{resolution.reason}",),
                ),
            )
        return (option,)

    for component_order in candidate_multistep_orders(component_ids):
        component_options = try_decompose_multistep_order(
            option=option,
            component_ids=component_order,
            client=client,
            ignored_common_compounds=ignored_common_compounds,
        )
        if component_options:
            return component_options

    return fallback_multistep_option(option)


def build_reaction_options(
    compound_id: str,
    client: KeggRestClient,
    ignored_common_compounds: set[str],
) -> List[ReactionOption]:
    """枚举所有能够以左→右或右→左方向生成指定化合物的候选反应。"""

    options: List[ReactionOption] = []
    reaction_ids = client.get_reaction_ids_for_compound(compound_id)
    client.prefetch_reactions(reaction_ids)

    for reaction_id in reaction_ids:
        reaction = client.try_get_reaction(reaction_id)
        if reaction is None:
            continue
        for direction, product_compounds in (
            ("left_to_right", stoichiometry_compound_ids(reaction.right_stoichiometry)),
            ("right_to_left", stoichiometry_compound_ids(reaction.left_stoichiometry)),
        ):
            if compound_id not in product_compounds:
                continue
            options.append(
                make_reaction_option(
                    compound_id=compound_id,
                    reaction=reaction,
                    direction=direction,
                    ignored_common_compounds=ignored_common_compounds,
                )
            )

    return options


def filter_reaction_options(
    options: Sequence[ReactionOption],
    allowed_reaction_ids: set[str] | None,
) -> List[ReactionOption]:
    """按可选的 module 反应白名单过滤候选；``None`` 表示不过滤。"""

    filtered: List[ReactionOption] = []
    for option in options:
        if not reaction_option_allowed_by_modules(option, allowed_reaction_ids):
            continue
        filtered.append(option)
    return filtered


def validate_electron_avoidance_mode(mode: str) -> str:
    """校验电子风险规避模式，避免静默使用未知策略。"""

    normalized = str(mode or DEFAULT_ELECTRON_AVOIDANCE_MODE).strip()
    if normalized not in ELECTRON_AVOIDANCE_MODES:
        raise ValueError(
            "electron_avoidance_mode must be one of: "
            + ", ".join(sorted(ELECTRON_AVOIDANCE_MODES))
            + f". Got: {mode}"
        )
    return normalized


def apply_electron_avoidance_filter(
    options: Sequence[ReactionOption],
    electron_avoidance_mode: str,
) -> List[ReactionOption]:
    """在严格模式下优先删除高电子风险候选。

    如果当前化合物的所有候选均为高风险，则全部保留，防止单个节点因电子
    风险直接断路；全局 ``strict_with_fallback`` 兜底逻辑在搜索调度层处理。
    """

    option_list = list(options)
    if electron_avoidance_mode not in {"strict", "strict_with_fallback"}:
        return option_list
    lower_risk = [
        option for option in option_list
        if option.electron_requirement.risk_score < HIGH_ELECTRON_RISK_SCORE
    ]
    if not lower_risk:
        return option_list
    high_risk = [
        option for option in option_list
        if option.electron_requirement.risk_score >= HIGH_ELECTRON_RISK_SCORE
    ]
    if not high_risk:
        return lower_risk
    # 保留一个高风险备用分支，防止低风险候选在更上游落入 incomplete
    # 死路后无法回溯到计量完整的显式电子反应。
    backup = min(
        high_risk,
        key=lambda option: (
            option.electron_requirement.risk_score,
            option.reaction.reaction_id,
            option.direction,
        ),
    )
    return [*lower_risk, backup]


# ---------------------------------------------------------------------------
# 底盘可达集合、内源方向索引和目标 ID 对齐
# ---------------------------------------------------------------------------

def load_reachable_compounds(path: Path) -> Tuple[set[str], pd.DataFrame]:
    """读取底盘可达性分析生成的 CSV，并返回合法 KEGG ID 集合。"""

    if not path.exists():
        raise FileNotFoundError(f"Reachable KEGG file not found: {path}")

    df = pd.read_csv(path)
    if "kegg_id" not in df.columns:
        raise ValueError(f"Missing kegg_id column in {path}")

    df["kegg_id"] = df["kegg_id"].astype(str).map(normalize_kegg_id)
    reachable = {kid for kid in df["kegg_id"] if KEGG_ID_PATTERN.match(kid)}
    return reachable, df


def wait_for_reachable_compounds_file(path: Path, timeout_seconds: int = REACHABLE_FILE_WAIT_SECONDS) -> None:
    """在底盘分析仍运行时等待其原子地产生可达化合物文件。

    running marker 不存在意味着上游没有运行或已经失败，此时立即报错；marker
    存在时才轮询等待，从而避免通路搜索读到尚未写完的中间文件。
    """

    if path.exists():
        return

    marker_path = path.parent / CHASSIS_ANALYSIS_RUNNING_MARKER
    startup_deadline = time.monotonic() + REACHABLE_FILE_STARTUP_GRACE_SECONDS
    while not marker_path.exists() and time.monotonic() < startup_deadline:
        if path.exists():
            return
        time.sleep(REACHABLE_MARKER_POLL_SECONDS)

    if not marker_path.exists():
        raise FileNotFoundError(f"Reachable KEGG file not found: {path}")

    deadline = time.monotonic() + timeout_seconds
    while time.monotonic() < deadline:
        if path.exists():
            return
        if not marker_path.exists():
            break
        time.sleep(REACHABLE_FILE_POLL_SECONDS)

    if path.exists():
        return
    raise FileNotFoundError(
        f"Reachable KEGG file not found after waiting for chassis analysis: {path}"
    )


def option_is_endogenous(
    option: ReactionOption,
    endogenous_reactions: set[str],
    endogenous_direction_index: EndogenousDirectionIndex | None = None,
) -> bool:
    """判断当前“反应 + 产物方向”是否属于底盘内源能力。

    首选方向索引：只有 GEM 上下界允许且该方向能够生成当前产物才算内源。
    旧调用方未提供方向索引时才退回仅按 reaction ID 判断。这里不会解析 GPR
    表达式，因此“内源”表示模型注释和边界支持，不等于已验证具体基因表达。
    """

    if endogenous_direction_index is not None:
        return is_directionally_endogenous(
            option.reaction.reaction_id,
            option.produced_compound,
            endogenous_direction_index,
        )
    return option.reaction.reaction_id in endogenous_reactions


# ---------------------------------------------------------------------------
# KEGG module 上下文构建
# ---------------------------------------------------------------------------

def resolve_target_module_context(
    target_compound: str,
    client: KeggRestClient,
    ignored_common_compounds: set[str],
) -> Tuple[Tuple[str, ...], Tuple[str, ...], Tuple[str, ...]]:
    """筛出真正覆盖目标生成反应的 KEGG module。

    化合物条目列出的 module 不一定都负责合成该化合物，因此先取直接生成目标
    的反应做交集；若 module 名称明确写成 ``... => 目标``，再优先保留这类
    生物合成 module。
    """

    compound_record = client.get_compound_record(target_compound)
    target_compound_module_ids = compound_record.module_ids

    direct_options = build_reaction_options(target_compound, client, ignored_common_compounds)
    direct_options = filter_reaction_options(direct_options, allowed_reaction_ids=None)
    direct_producing_reaction_ids = dedupe_keep_order(
        option.reaction.reaction_id for option in direct_options
    )
    direct_reaction_id_set = set(direct_producing_reaction_ids)

    client.prefetch_modules(target_compound_module_ids)
    matched_target_module_ids = dedupe_keep_order(
        module_id
        for module_id in target_compound_module_ids
        if direct_reaction_id_set & set(client.get_module(module_id).reaction_ids)
        or any(
            reaction_option_allowed_by_modules(option, set(client.get_module(module_id).reaction_ids))
            for option in direct_options
        )
    )
    synthesis_named_module_ids = dedupe_keep_order(
        module_id
        for module_id in matched_target_module_ids
        if module_name_targets_compound(client.get_module(module_id).name, compound_record.name)
    )
    if synthesis_named_module_ids:
        matched_target_module_ids = synthesis_named_module_ids
    return (
        target_compound_module_ids,
        matched_target_module_ids,
        direct_producing_reaction_ids,
    )


def find_upstream_modules_for_boundary_compound(
    compound_id: str,
    client: KeggRestClient,
    ignored_common_compounds: set[str],
) -> Tuple[str, ...]:
    """查找能够生成某个 module 边界化合物的上游 module。"""

    module_ids: List[str] = []
    options = build_reaction_options(compound_id, client, ignored_common_compounds)
    for option in options:
        for module_id in option.reaction.module_ids:
            module = client.get_module(module_id)
            if compound_id in module.produced_compound_ids or reaction_option_allowed_by_modules(
                option,
                set(module.reaction_ids),
            ):
                module_ids.append(module_id)
    return dedupe_keep_order(module_ids)


def find_unresolved_module_start_compounds(
    module: ModuleRecord,
    seed_compound_ids: Sequence[str],
    reachable_compounds: set[str],
    ignored_common_compounds: set[str],
) -> Tuple[str, ...]:
    """在单个 module 内反向追踪，找出仍未连接到底盘的入口化合物。"""

    pending = list(
        compound_id
        for compound_id in dedupe_keep_order(seed_compound_ids)
        if compound_id not in reachable_compounds
        and compound_id not in ignored_common_compounds
    )
    visited = set(pending)

    while pending:
        current_compound = pending.pop(0)
        producing_edges = [
            edge
            for edge in module.reaction_edges
            if current_compound in edge.produced_compound_ids
        ]
        for edge in producing_edges:
            for precursor_id in edge.consumed_compound_ids:
                if precursor_id in ignored_common_compounds:
                    continue
                if precursor_id in reachable_compounds:
                    continue
                if precursor_id in visited:
                    continue
                visited.add(precursor_id)
                pending.append(precursor_id)

    module_start_compound_ids = set(module.start_compound_ids)
    return dedupe_keep_order(
        compound_id
        for compound_id in visited
        if compound_id in module_start_compound_ids
        and compound_id not in reachable_compounds
        and compound_id not in ignored_common_compounds
    )


def build_module_chain_context(
    initial_module_ids: Sequence[str],
    seed_compound_ids: Sequence[str],
    reachable_compounds: set[str],
    client: KeggRestClient,
    ignored_common_compounds: set[str],
    max_module_chain_depth: int,
) -> ModuleChainContext:
    """从目标 module 向上游递归拼接 module，并汇总允许反应白名单。"""

    expanded_module_ids: List[str] = list(dedupe_keep_order(initial_module_ids))
    boundary_compound_ids: List[str] = []
    seen_modules = set(expanded_module_ids)
    frontier_module_seeds = tuple(
        (module_id, tuple(dedupe_keep_order(seed_compound_ids)))
        for module_id in expanded_module_ids
    )

    # 每一轮以当前 module 的未解决入口为边界，再寻找能够生成边界的上游 module。
    for _ in range(max(0, max_module_chain_depth)):
        frontier_module_ids = tuple(module_id for module_id, _ in frontier_module_seeds)
        client.prefetch_modules(frontier_module_ids)
        frontier_boundary_compounds = dedupe_keep_order(
            compound_id
            for module_id, module_seed_compound_ids in frontier_module_seeds
            for compound_id in find_unresolved_module_start_compounds(
                module=client.get_module(module_id),
                seed_compound_ids=module_seed_compound_ids,
                reachable_compounds=reachable_compounds,
                ignored_common_compounds=ignored_common_compounds,
            )
        )
        if not frontier_boundary_compounds:
            break

        boundary_compound_ids.extend(frontier_boundary_compounds)
        next_frontier_module_seeds: List[Tuple[str, Tuple[str, ...]]] = []
        for compound_id in frontier_boundary_compounds:
            for module_id in find_upstream_modules_for_boundary_compound(
                compound_id,
                client,
                ignored_common_compounds,
            ):
                if module_id in seen_modules:
                    continue
                module = client.get_module(module_id)
                if not module.reaction_edges:
                    continue
                seen_modules.add(module_id)
                expanded_module_ids.append(module_id)
                next_frontier_module_seeds.append((module_id, (compound_id,)))

        frontier_module_seeds = tuple(next_frontier_module_seeds)
        if not frontier_module_seeds:
            break

    allowed_reaction_ids = dedupe_keep_order(
        reaction_id
        for module_id in expanded_module_ids
        for reaction_id in client.get_module(module_id).reaction_ids
    )
    return ModuleChainContext(
        expanded_module_ids=tuple(expanded_module_ids),
        boundary_compound_ids=dedupe_keep_order(boundary_compound_ids),
        allowed_reaction_ids=allowed_reaction_ids,
    )


def normalize_text_for_matching(text: str) -> str:
    """把名称规范为小写字母数字词串，供宽松包含匹配使用。"""

    return re.sub(r"[^a-z0-9]+", " ", text.lower()).strip()


def module_name_targets_compound(module_name: str, compound_name: str) -> bool:
    """判断 module 名称箭头右侧是否明确指向目标化合物。"""

    if "=>" not in module_name:
        return False
    _, _, right_side = module_name.partition("=>")
    normalized_compound = normalize_text_for_matching(compound_name)
    normalized_right = normalize_text_for_matching(right_side)
    return bool(normalized_compound) and normalized_compound in normalized_right


def choose_expand_compound(
    unresolved: frozenset[str],
    client: KeggRestClient,
) -> str:
    """选择候选反应最少的未解决化合物，采用 fail-first 策略缩小分支。"""

    return min(
        unresolved,
        key=lambda compound_id: (len(client.get_reaction_ids_for_compound(compound_id)), compound_id),
    )


# ---------------------------------------------------------------------------
# 候选排序、搜索状态与剪枝
# ---------------------------------------------------------------------------

def rank_reaction_options(
    options: Sequence[ReactionOption],
    reachable_compounds: set[str],
    endogenous_reactions: set[str],
    client: KeggRestClient,
    ignored_common_compounds: set[str],
    target_pathway_ids: set[str],
    max_major_precursors: int,
    allowed_reaction_ids: set[str] | None = None,
    electron_avoidance_mode: str = "prefer",
    endogenous_direction_index: EndogenousDirectionIndex | None = None,
) -> List[ReactionOption]:
    """按生物学启发式对生成当前化合物的反应候选做字典序排序。

    排序综合考虑反应方向、SAM/NADPH/CoA 负担、电子风险、前体与底盘
    可达集合的距离、同通路证据和酶注释完整性。这里仅改变探索顺序，不会把
    ``disfavored`` 候选直接判为不可行；最终可行性仍交给下游 GEM 验证。
    """

    electron_avoidance_mode = validate_electron_avoidance_mode(electron_avoidance_mode)
    one_step_cache: Dict[str, int] = {}

    def estimate_distance_to_reachable(compound_id: str) -> int:
        """用至多一层前瞻估计前体距可达集合的离散距离（0～3）。"""

        if compound_id in one_step_cache:
            return one_step_cache[compound_id]
        if compound_id in reachable_compounds:
            one_step_cache[compound_id] = 0
            return 0

        best = 3
        compound_options = build_reaction_options(compound_id, client, ignored_common_compounds)
        compound_options = filter_reaction_options(compound_options, allowed_reaction_ids)
        for option in compound_options[:8]:
            precursor_set = set(option.precursor_compounds)
            if len(precursor_set) > max_major_precursors:
                continue
            if not precursor_set:
                best = min(best, 1)
                continue
            if all(precursor in reachable_compounds for precursor in precursor_set):
                best = min(best, 1)
                continue
            if any(precursor in reachable_compounds for precursor in precursor_set):
                best = min(best, 2)

        one_step_cache[compound_id] = best
        return best

    # 先控制主要前体分支数，再按电子规避模式做局部过滤，限制组合爆炸。
    filtered: List[ReactionOption] = []

    for option in options:
        precursor_set = set(option.precursor_compounds)
        if len(precursor_set) > max_major_precursors:
            continue
        filtered.append(option)
    filtered = apply_electron_avoidance_filter(filtered, electron_avoidance_mode)

    scored: List[Tuple[Tuple[Any, ...], ReactionOption]] = []
    for option in filtered:
        precursor_set = set(option.precursor_compounds)

        reachable_count = sum(1 for compound_id in precursor_set if compound_id in reachable_compounds)
        unresolved_count = sum(1 for compound_id in precursor_set if compound_id not in reachable_compounds)
        distance_score = sum(
            estimate_distance_to_reachable(compound_id)
            for compound_id in precursor_set
            if compound_id not in reachable_compounds
        )
        thermo_penalty = {"favored": 0, "neutral": 1, "disfavored": 2}.get(
            option.screening.thermo_direction,
            1,
        )
        endogenous_penalty = 0 if option_is_endogenous(
            option,
            endogenous_reactions,
            endogenous_direction_index,
        ) else 1
        shared_pathway_penalty = 0 if target_pathway_ids & set(option.reaction.pathway_ids) else 1
        enzyme_penalty = 0 if option.reaction.enzyme_ecs or option.reaction.ko_ids else 1

        reachable_penalty = 0 if reachable_count > 0 else 1
        electron_risk_penalty = (
            0
            if electron_avoidance_mode == "off"
            else option.electron_requirement.risk_score
        )
        # 元组按字段依次比较，因此前面的方向和辅因子指标优先级更高。
        score = (
            thermo_penalty,
            option.screening.sam_burden,
            option.screening.nadph_burden,
            option.screening.coa_burden,
            electron_risk_penalty,
            reachable_penalty,
            distance_score,
            unresolved_count,
            shared_pathway_penalty,
            enzyme_penalty,
            endogenous_penalty,
            len(precursor_set),
            option.reaction.reaction_id,
            option.direction,
        )
        scored.append((score, option))

    scored.sort(key=lambda item: item[0])
    return [option for _, option in scored]


def canonical_state(unresolved: Iterable[str]) -> frozenset[str]:
    """将未解决前体集合规范为可哈希、与遍历顺序无关的状态键。"""

    return frozenset(sorted(unresolved))


def plan_route_signature(plan_steps: Sequence[PlanStep]) -> StateRouteSignature:
    """将路线历史压缩为用于去重的“产物—反应—方向”签名。"""

    return tuple(
        (
            step.option.produced_compound,
            step.option.reaction.reaction_id,
            step.option.direction,
        )
        for step in plan_steps
    )


def state_route_sort_key(
    record: StateRouteRecord,
) -> Tuple[int, int, int, int, StateRouteSignature]:
    """定义同一状态下历史路线的保留优先级。"""

    max_electron_risk_score, heterologous_steps, total_steps, route_signature = record
    return (
        max_electron_risk_score,
        heterologous_steps,
        total_steps,
        len(route_signature),
        route_signature,
    )


def remember_state_route(
    routes_by_state: Dict[frozenset[str], List[StateRouteRecord]],
    state: frozenset[str],
    heterologous_steps: int,
    total_steps: int,
    route_signature: StateRouteSignature,
    max_routes_per_state: int,
    max_electron_risk_score: int = 0,
) -> bool:
    """为同一未解决状态保留成本最低且历史不同的前 K 条路线。

    不能只为每个状态保留一条路线：两条路线即使得到相同未解决前体，也可能
    已使用不同反应，后续受到不同的环路约束。返回值表示新路线是否进入 Top-K。
    """

    route_limit = max(1, max_routes_per_state)
    for record in routes_by_state.get(state, []):
        if record[3] == route_signature and record[:3] <= (
            max_electron_risk_score,
            heterologous_steps,
            total_steps,
        ):
            return False

    records = [
        record
        for record in routes_by_state.get(state, [])
        if record[3] != route_signature
    ]
    records.append(
        (
            max_electron_risk_score,
            heterologous_steps,
            total_steps,
            route_signature,
        )
    )
    records.sort(key=state_route_sort_key)
    kept = records[:route_limit]
    routes_by_state[state] = kept
    return any(
        record[0] == max_electron_risk_score
        and record[1] == heterologous_steps
        and record[2] == total_steps
        and record[3] == route_signature
        for record in kept
    )


def is_state_route_active(
    routes_by_state: Dict[frozenset[str], List[StateRouteRecord]],
    state: frozenset[str],
    heterologous_steps: int,
    total_steps: int,
    route_signature: StateRouteSignature,
    max_electron_risk_score: int = 0,
) -> bool:
    """检查从堆中弹出的路线是否仍在该状态的 Top-K 活跃集合中。"""

    return any(
        record[0] == max_electron_risk_score
        and record[1] == heterologous_steps
        and record[2] == total_steps
        and record[3] == route_signature
        for record in routes_by_state.get(state, [])
    )


def split_precursors_by_reachability(
    precursors: Sequence[str],
    reachable_compounds: set[str],
) -> Tuple[Tuple[str, ...], Tuple[str, ...]]:
    """把主要前体拆为底盘已可达和仍需继续反向搜索的两组。"""

    reachable_precursors = tuple(compound_id for compound_id in precursors if compound_id in reachable_compounds)
    unresolved_precursors = tuple(compound_id for compound_id in precursors if compound_id not in reachable_compounds)
    return reachable_precursors, unresolved_precursors


def expand_unresolved_state(
    current_compound: str,
    unresolved: frozenset[str],
    option: ReactionOption,
    reachable_compounds: set[str],
) -> Tuple[frozenset[str], Tuple[str, ...], Tuple[str, ...]]:
    """用一个反应替换当前目标，得到下一未解决状态及其可达性证据。"""

    reachable_precursors, unresolved_precursors = split_precursors_by_reachability(
        option.precursor_compounds,
        reachable_compounds,
    )
    next_unresolved = set(unresolved)
    next_unresolved.discard(current_compound)
    next_unresolved.update(unresolved_precursors)
    return canonical_state(next_unresolved), reachable_precursors, unresolved_precursors


def make_plan_step(
    option: ReactionOption,
    reachable_precursors: Tuple[str, ...],
    is_endogenous: bool,
) -> PlanStep:
    """把反应选项和本次扩展结果冻结为可输出的路线步骤。"""

    return PlanStep(
        option=option,
        reachable_precursors=reachable_precursors,
        is_endogenous=is_endogenous,
        source_reaction_ids=(option.reaction.reaction_id,),
        resolution_action=(
            "explicit_multistep_component"
            if any(
                hit.startswith("decomposed_from:")
                for hit in option.screening.rule_hits
            )
            else "none"
        ),
        resolution_evidence=tuple(
            hit
            for hit in option.screening.rule_hits
            if hit.startswith((
                "decomposed_from:",
                "decomposition_index:",
                "reaction_resolution:",
                "multistep_decomposition_fallback:",
            ))
        ),
    )


def validate_reaction_resolution_mode(mode: str) -> str:
    normalized = str(mode or "strict").strip().lower()
    if normalized not in REACTION_RESOLUTION_MODES:
        raise ValueError(
            "reaction_resolution_mode must be one of: "
            + ", ".join(sorted(REACTION_RESOLUTION_MODES))
            + f". Got: {mode}"
        )
    return normalized


def _net_stoichiometry_for_option(
    option: ReactionOption,
    ignored_common_compounds: set[str],
) -> Dict[str, float]:
    """计算一个定向反应的非通用化合物净计量，生成物为正。"""

    consumed, produced = directional_stoichiometry(option.reaction, option.direction)
    net: Dict[str, float] = {}
    for compound_id, amount in consumed:
        if compound_id not in ignored_common_compounds:
            net[compound_id] = net.get(compound_id, 0.0) - float(amount)
    for compound_id, amount in produced:
        if compound_id not in ignored_common_compounds:
            net[compound_id] = net.get(compound_id, 0.0) + float(amount)
    return {
        compound_id: amount
        for compound_id, amount in net.items()
        if abs(amount) > 1e-9
    }


def _net_signature(net: Dict[str, float]) -> Tuple[Tuple[str, float], ...]:
    return tuple(sorted(
        (compound_id, round(amount, 9))
        for compound_id, amount in net.items()
        if abs(amount) > 1e-9
    ))


def _combined_net_signature(
    steps: Sequence[PlanStep],
    ignored_common_compounds: set[str],
) -> Tuple[Tuple[str, float], ...]:
    net: Dict[str, float] = {}
    for step in steps:
        for compound_id, amount in _net_stoichiometry_for_option(
            step.option,
            ignored_common_compounds,
        ).items():
            net[compound_id] = net.get(compound_id, 0.0) + amount
    return _net_signature(net)


def _adjacent_pair_target(
    first: PlanStep,
    second: PlanStep,
) -> str | None:
    """返回相邻反向搜索步骤中更下游的产物；无唯一衔接时返回 None。"""

    first_uses_second = second.option.produced_compound in first.option.precursor_compounds
    second_uses_first = first.option.produced_compound in second.option.precursor_compounds
    if first_uses_second == second_uses_first:
        return None
    return (
        first.option.produced_compound
        if first_uses_second
        else second.option.produced_compound
    )


def _find_unique_complete_pair_replacement(
    first: PlanStep,
    second: PlanStep,
    client: KeggRestClient,
    ignored_common_compounds: set[str],
) -> ReactionOption | None:
    """用净计量匹配为相邻不完整步骤寻找唯一完整 KEGG 反应。"""

    target_compound = _adjacent_pair_target(first, second)
    if not target_compound:
        return None
    pair_signature = _combined_net_signature(
        (first, second),
        ignored_common_compounds,
    )
    source_ids = {
        first.option.reaction.reaction_id,
        second.option.reaction.reaction_id,
    }
    candidates: Dict[Tuple[str, str], ReactionOption] = {}
    for option in build_reaction_options(
        target_compound,
        client,
        ignored_common_compounds,
    ):
        reaction = option.reaction
        if reaction.reaction_id in source_ids:
            continue
        resolution = classify_reaction_resolution(reaction)
        if resolution.is_incomplete or resolution.hard_blocker:
            continue
        if not reaction.enzyme_ecs and not reaction.ko_ids:
            continue
        if _net_signature(
            _net_stoichiometry_for_option(option, ignored_common_compounds)
        ) != pair_signature:
            continue
        candidates[(reaction.reaction_id, option.direction)] = option
    if len(candidates) != 1:
        return None
    return next(iter(candidates.values()))


def normalize_solution_reactions(
    solution: Solution,
    client: KeggRestClient,
    reachable_compounds: set[str],
    endogenous_reactions: set[str],
    ignored_common_compounds: set[str],
    endogenous_direction_index: EndogenousDirectionIndex | None,
) -> Solution:
    """合并可唯一解析的相邻不完整步骤，并标记仍不可实现的反应。"""

    steps = list(solution.steps)
    events: List[str] = []
    changed = True
    while changed:
        changed = False
        for index in range(len(steps) - 1):
            first, second = steps[index], steps[index + 1]
            if not (
                classify_reaction_resolution(first.option.reaction).is_incomplete
                or classify_reaction_resolution(second.option.reaction).is_incomplete
            ):
                continue
            replacement = _find_unique_complete_pair_replacement(
                first,
                second,
                client,
                ignored_common_compounds,
            )
            if replacement is None:
                continue
            source_ids = dedupe_keep_order((
                *first.source_reaction_ids,
                *second.source_reaction_ids,
            ))
            merge_event = (
                f"merged_from:{'+'.join(source_ids)}"
                f"->{replacement.reaction.reaction_id}"
            )
            replacement = add_screening_hits(
                replacement,
                (merge_event, "reaction_resolution:net_stoichiometry_match"),
            )
            reachable_precursors, _ = split_precursors_by_reachability(
                replacement.precursor_compounds,
                reachable_compounds,
            )
            replacement_step = PlanStep(
                option=replacement,
                reachable_precursors=reachable_precursors,
                is_endogenous=option_is_endogenous(
                    replacement,
                    endogenous_reactions,
                    endogenous_direction_index,
                ),
                source_reaction_ids=source_ids,
                resolution_action="merged_complete_reaction",
                resolution_evidence=(
                    merge_event,
                    "net_stoichiometry_match",
                    "unique_complete_kegg_reaction",
                ),
            )
            steps[index : index + 2] = [replacement_step]
            events.append(merge_event)
            changed = True
            break

    blocking_ids: List[str] = []
    blocking_reasons: List[str] = []
    for step in steps:
        resolution = classify_reaction_resolution(step.option.reaction)
        if not resolution.hard_blocker:
            continue
        blocking_ids.append(step.option.reaction.reaction_id)
        blocking_reasons.append(
            f"{step.option.reaction.reaction_id}:{resolution.reason}"
        )
    blocking_ids_tuple = dedupe_keep_order(blocking_ids)
    status = "blocked" if blocking_ids_tuple else (
        "resolved_with_merges" if events else "resolved"
    )
    return Solution(
        steps=tuple(steps),
        reaction_resolution_status=status,
        normalization_events=tuple(events),
        blocking_reaction_ids=blocking_ids_tuple,
        blocking_reasons=dedupe_keep_order(blocking_reasons),
    )


# ---------------------------------------------------------------------------
# 单轮反向搜索与分层兜底调度
# ---------------------------------------------------------------------------

def search_gap_solutions_once(
    target_compound: str,
    reachable_compounds: set[str],
    endogenous_reactions: set[str],
    client: KeggRestClient,
    max_total_steps: int,
    max_new_enzymes: int,
    max_solutions: int,
    max_reactions_per_compound: int,
    max_major_precursors: int,
    max_routes_per_state: int,
    ignored_common_compounds: set[str],
    allowed_reaction_ids: set[str] | None,
    electron_avoidance_mode: str,
    reaction_resolution_mode: str = "strict",
    endogenous_direction_index: EndogenousDirectionIndex | None = None,
) -> SearchRoundResult:
    """在固定 module 白名单、电子策略和资源上限下执行一次反向最佳优先搜索。

    状态是尚未接到底盘可达集合的主要前体集合。每次选择其中候选最少的一个
    化合物，用能够生成它的反应替换为更上游前体；当集合为空时得到完整路线。
    堆优先异源步骤较少、预计剩余工作较小的节点。
    """

    electron_avoidance_mode = validate_electron_avoidance_mode(electron_avoidance_mode)
    reaction_resolution_mode = validate_reaction_resolution_mode(
        reaction_resolution_mode
    )
    target_entry_pathways = set(client.get_compound_record(target_compound).pathway_ids)

    # 目标已经可达时返回零步路线，使调用方仍能得到结构一致的成功结果。
    if target_compound in reachable_compounds:
        return (
            Solution(
                steps=tuple(),
            ),
        )

    # 优先队列保存待扩展路线；计数器确保所有成本相同时仍能稳定排序。
    heap: List[SearchQueueItem] = []
    push_counter = itertools.count()
    start_state = canonical_state({target_compound})
    heapq.heappush(
        heap,
        SearchQueueItem(
            heterologous_steps=0,
            estimated_total_work=1,
            total_steps=0,
            unresolved_count=1,
            sequence=next(push_counter),
            unresolved=start_state,
            plan_steps=tuple(),
        ),
    )

    # 同一未解决集合保留多条不同历史，兼顾剪枝强度与环路约束下的路线多样性。
    routes_by_state: Dict[frozenset[str], List[StateRouteRecord]] = {
        start_state: [(0, 0, 0, tuple())],
    }
    solutions: List[Solution] = []
    rejected_solutions: List[Solution] = []
    rejected_signatures: set[StateRouteSignature] = set()

    while heap and len(solutions) < max_solutions:
        item = heapq.heappop(heap)
        heterologous_steps = item.heterologous_steps
        total_steps = item.total_steps
        unresolved = item.unresolved
        plan_steps = item.plan_steps
        route_signature = plan_route_signature(plan_steps)
        max_electron_risk_score = (
            0
            if electron_avoidance_mode == "off"
            else max(
                (
                    step.option.electron_requirement.risk_score
                    for step in plan_steps
                ),
                default=0,
            )
        )
        used_reaction_ids = {
            step.option.reaction.reaction_id
            for step in plan_steps
        }
        expanded_compounds = {
            step.option.produced_compound
            for step in plan_steps
        }

        # Top-K 集合可能在节点入堆后被更优路线更新，弹出时需丢弃过期节点。
        if unresolved and not is_state_route_active(
            routes_by_state=routes_by_state,
            state=unresolved,
            heterologous_steps=heterologous_steps,
            total_steps=total_steps,
            route_signature=route_signature,
            max_electron_risk_score=max_electron_risk_score,
        ):
            continue

        # 所有主要前体均已落入可达集合，当前路线即为一个完整解。
        if not unresolved:
            normalized_solution = normalize_solution_reactions(
                solution=Solution(steps=plan_steps),
                client=client,
                reachable_compounds=reachable_compounds,
                endogenous_reactions=endogenous_reactions,
                ignored_common_compounds=ignored_common_compounds,
                endogenous_direction_index=endogenous_direction_index,
            )
            if not normalized_solution.reaction_ready:
                if (
                    route_signature not in rejected_signatures
                    and len(rejected_solutions) < max(10, max_solutions * 5)
                ):
                    rejected_signatures.add(route_signature)
                    rejected_solutions.append(normalized_solution)
                # strict 和 audit 都不把硬阻断路线交给推荐/GEM；audit 的区别是
                # 输出保留更完整的拒绝证据，避免问题路线占满 max_solutions。
                continue
            solutions.append(normalized_solution)
            continue

        if total_steps >= max_total_steps:
            continue

        # 采用 fail-first 策略选节点，并只扩展排序最前的有限个候选反应。
        current_compound = choose_expand_compound(unresolved, client)
        all_options = build_reaction_options(current_compound, client, ignored_common_compounds)
        all_options = filter_reaction_options(all_options, allowed_reaction_ids)
        ranked_options = rank_reaction_options(
            options=all_options,
            reachable_compounds=reachable_compounds,
            endogenous_reactions=endogenous_reactions,
            client=client,
            ignored_common_compounds=ignored_common_compounds,
            target_pathway_ids=target_entry_pathways,
            max_major_precursors=max_major_precursors,
            allowed_reaction_ids=allowed_reaction_ids,
            electron_avoidance_mode=electron_avoidance_mode,
            endogenous_direction_index=endogenous_direction_index,
        )[:max_reactions_per_compound]

        for option in ranked_options:
            # KEGG 中部分反应是多步汇总条目；能可靠拆解时按真实组件逐步计数。
            component_options = decompose_multistep_option(
                option=option,
                client=client,
                ignored_common_compounds=ignored_common_compounds,
            )
            if not component_options:
                continue

            # 禁止自依赖、组件内重复反应和路线历史中的反应复用，避免显式环路。
            if any(component_option.produced_compound in component_option.precursor_compounds for component_option in component_options):
                continue
            component_reaction_ids = [
                component_option.reaction.reaction_id
                for component_option in component_options
            ]
            if len(component_reaction_ids) != len(set(component_reaction_ids)):
                continue
            if used_reaction_ids.intersection(component_reaction_ids):
                continue

            next_state = unresolved
            next_plan_steps: List[PlanStep] = []
            component_heterologous_steps = 0
            expansion_blocked = False
            component_expanded_compounds = set(expanded_compounds)
            for component_option in component_options:
                previous_state = next_state
                next_state, reachable_precursors, unresolved_precursors = expand_unresolved_state(
                    current_compound=component_option.produced_compound,
                    unresolved=next_state,
                    option=component_option,
                    reachable_compounds=reachable_compounds,
                )
                # 状态不变表示该步没有解决任何目标；重新引入已扩展产物会形成回路。
                if next_state == previous_state:
                    expansion_blocked = True
                    break
                if component_expanded_compounds.intersection(unresolved_precursors):
                    expansion_blocked = True
                    break

                is_endogenous = option_is_endogenous(
                    component_option,
                    endogenous_reactions,
                    endogenous_direction_index,
                )
                component_heterologous_steps += 0 if is_endogenous else 1
                next_plan_steps.append(
                    make_plan_step(
                        reachable_precursors=reachable_precursors,
                        option=component_option,
                        is_endogenous=is_endogenous,
                    )
                )
                component_expanded_compounds.add(component_option.produced_compound)

            if expansion_blocked:
                continue
            if next_state == unresolved:
                continue

            next_heterologous_steps = heterologous_steps + component_heterologous_steps
            next_total_steps = total_steps + len(next_plan_steps)

            # 每个未解决化合物至少还需一步，据此做乐观下界剪枝。
            if next_heterologous_steps > max_new_enzymes:
                continue
            if next_total_steps + len(next_state) > max_total_steps:
                continue

            next_plan = plan_steps + tuple(next_plan_steps)
            next_route_signature = plan_route_signature(next_plan)
            next_max_electron_risk_score = (
                0
                if electron_avoidance_mode == "off"
                else max(
                    max_electron_risk_score,
                    *(
                        step.option.electron_requirement.risk_score
                        for step in next_plan_steps
                    ),
                )
            )
            # 只有进入该状态 Top-K 的历史才继续入堆，控制状态组合爆炸。
            if next_state and not remember_state_route(
                routes_by_state=routes_by_state,
                state=next_state,
                heterologous_steps=next_heterologous_steps,
                total_steps=next_total_steps,
                route_signature=next_route_signature,
                max_routes_per_state=max_routes_per_state,
                max_electron_risk_score=next_max_electron_risk_score,
            ):
                continue

            heapq.heappush(
                heap,
                SearchQueueItem(
                    heterologous_steps=next_heterologous_steps,
                    estimated_total_work=next_total_steps + len(next_state),
                    total_steps=next_total_steps,
                    unresolved_count=len(next_state),
                    sequence=next(push_counter),
                    unresolved=next_state,
                    plan_steps=next_plan,
                ),
            )

    return SearchRoundResult(
        solutions=tuple(solutions),
        rejected_solutions=tuple(rejected_solutions),
    )


def run_search_gap_analysis(
    target_compound: str,
    reachable_compounds: set[str],
    endogenous_reactions: set[str],
    client: KeggRestClient,
    max_total_steps: int,
    max_new_enzymes: int,
    max_solutions: int,
    max_reactions_per_compound: int,
    max_major_precursors: int,
    ignored_common_compounds: set[str],
    module_filter_mode: str = DEFAULT_MODULE_FILTER_MODE,
    max_routes_per_state: int = DEFAULT_MAX_ROUTES_PER_STATE,
    max_module_chain_depth: int = DEFAULT_MAX_MODULE_CHAIN_DEPTH,
    electron_avoidance_mode: str = DEFAULT_ELECTRON_AVOIDANCE_MODE,
    reaction_resolution_mode: str = "strict",
    endogenous_direction_index: EndogenousDirectionIndex | None = None,
) -> SearchExecutionResult:
    """构建 module 上下文并调度 module 与电子风险两层搜索兜底。

    ``module_filter_mode`` 控制是否先在 KEGG module 链内搜索；``prefer`` 在
    module 内无解时回退到完整反应空间。``strict_with_fallback`` 则先规避高
    电子风险候选，无解后放宽为 prefer，并同时扩大搜索资源上限。
    """

    electron_avoidance_mode = validate_electron_avoidance_mode(electron_avoidance_mode)
    reaction_resolution_mode = validate_reaction_resolution_mode(
        reaction_resolution_mode
    )
    (
        target_compound_module_ids,
        matched_target_module_ids,
        direct_producing_reaction_ids,
    ) = resolve_target_module_context(target_compound, client, ignored_common_compounds)

    # 关闭 module 筛选时不构造白名单；否则从目标 module 向上游扩展反应子图。
    if module_filter_mode == "off":
        module_chain_context = ModuleChainContext(
            expanded_module_ids=tuple(),
            boundary_compound_ids=tuple(),
            allowed_reaction_ids=tuple(),
        )
    else:
        module_chain_context = build_module_chain_context(
            initial_module_ids=matched_target_module_ids,
            seed_compound_ids=(target_compound,),
            reachable_compounds=reachable_compounds,
            client=client,
            ignored_common_compounds=ignored_common_compounds,
            max_module_chain_depth=max_module_chain_depth,
        )
    allowed_reaction_id_set = set(module_chain_context.allowed_reaction_ids)
    # 基础上限来自调用参数；电子兜底只会放宽，不会缩小用户给定上限。
    base_search_limits = {
        "max_total_steps": max_total_steps,
        "max_new_enzymes": max_new_enzymes,
        "max_reactions_per_compound": max_reactions_per_compound,
        "max_routes_per_state": max_routes_per_state,
    }
    expanded_fallback_limits = {
        "max_total_steps": max(max_total_steps, DEFAULT_ELECTRON_FALLBACK_MAX_TOTAL_STEPS),
        "max_new_enzymes": max(max_new_enzymes, DEFAULT_ELECTRON_FALLBACK_MAX_NEW_ENZYMES),
        "max_reactions_per_compound": max(
            max_reactions_per_compound,
            DEFAULT_ELECTRON_FALLBACK_MAX_REACTIONS_PER_COMPOUND,
        ),
        "max_routes_per_state": max(
            max_routes_per_state,
            DEFAULT_ELECTRON_FALLBACK_MAX_ROUTES_PER_STATE,
        ),
    }

    def run_once(
        allowed_reaction_ids: set[str] | None,
        active_electron_mode: str,
        search_limits: Dict[str, int],
    ) -> SearchRoundResult:
        """把本阶段白名单、电子策略和资源上限传给单轮搜索。"""

        return search_gap_solutions_once(
            target_compound=target_compound,
            reachable_compounds=reachable_compounds,
            endogenous_reactions=endogenous_reactions,
            client=client,
            max_total_steps=search_limits["max_total_steps"],
            max_new_enzymes=search_limits["max_new_enzymes"],
            max_solutions=max_solutions,
            max_reactions_per_compound=search_limits["max_reactions_per_compound"],
            max_major_precursors=max_major_precursors,
            max_routes_per_state=search_limits["max_routes_per_state"],
            ignored_common_compounds=ignored_common_compounds,
            allowed_reaction_ids=allowed_reaction_ids,
            electron_avoidance_mode=active_electron_mode,
            reaction_resolution_mode=reaction_resolution_mode,
            endogenous_direction_index=endogenous_direction_index,
        )

    def run_module_pipeline(
        active_electron_mode: str,
        search_limits: Dict[str, int],
    ) -> Tuple[SearchRoundResult, str, bool]:
        """执行 module 筛选策略，并返回解、实际模式及是否发生全空间回退。"""

        search_mode_used = "full_only"
        did_fallback = False
        round_result = SearchRoundResult(tuple(), tuple())
        if module_filter_mode == "off":
            round_result = run_once(None, active_electron_mode, search_limits)
        elif module_filter_mode == "strict":
            search_mode_used = "module_chained_strict"
            if module_chain_context.expanded_module_ids:
                round_result = run_once(allowed_reaction_id_set, active_electron_mode, search_limits)
        else:
            if module_chain_context.expanded_module_ids:
                search_mode_used = "module_chained"
                round_result = run_once(allowed_reaction_id_set, active_electron_mode, search_limits)

            if not round_result.solutions:
                did_fallback = bool(module_chain_context.expanded_module_ids)
                search_mode_used = "fallback_full" if did_fallback else "full_only"
                fallback_result = run_once(None, active_electron_mode, search_limits)
                round_result = SearchRoundResult(
                    solutions=fallback_result.solutions,
                    rejected_solutions=(
                        *round_result.rejected_solutions,
                        *fallback_result.rejected_solutions,
                    ),
                )
        return round_result, search_mode_used, did_fallback

    # strict_with_fallback 的第二阶段同时放宽电子偏好和搜索宽度/深度。
    electron_avoidance_fallback = False
    electron_fallback_parameters: Dict[str, int] = {}
    active_electron_mode = "strict" if electron_avoidance_mode == "strict_with_fallback" else electron_avoidance_mode
    round_result, search_mode_used, did_fallback = run_module_pipeline(active_electron_mode, base_search_limits)
    if electron_avoidance_mode == "strict_with_fallback" and not round_result.solutions:
        electron_avoidance_fallback = True
        electron_fallback_parameters = dict(expanded_fallback_limits)
        fallback_round, search_mode_used, did_fallback = run_module_pipeline("prefer", expanded_fallback_limits)
        round_result = SearchRoundResult(
            solutions=fallback_round.solutions,
            rejected_solutions=(
                *round_result.rejected_solutions,
                *fallback_round.rejected_solutions,
            ),
        )
        search_mode_used = f"electron_expanded_fallback_{search_mode_used}"

    return SearchExecutionResult(
        solutions=round_result.solutions,
        rejected_solutions=round_result.rejected_solutions,
        search_mode_used=search_mode_used,
        did_fallback=did_fallback,
        electron_avoidance_mode=electron_avoidance_mode,
        electron_avoidance_fallback=electron_avoidance_fallback,
        electron_fallback_parameters=electron_fallback_parameters,
        target_compound_module_ids=target_compound_module_ids,
        matched_target_module_ids=matched_target_module_ids,
        direct_producing_reaction_ids=direct_producing_reaction_ids,
        module_chain_context=module_chain_context,
        reaction_resolution_mode=reaction_resolution_mode,
    )


# ---------------------------------------------------------------------------
# 结果汇总与持久化
# ---------------------------------------------------------------------------

def compound_label(compound_id: str, client: KeggRestClient) -> str:
    """生成“KEGG ID (名称)”标签；名称获取失败时仅保留 ID。"""

    try:
        name = client.get_compound_name(compound_id)
    except Exception:
        name = ""
    return f"{compound_id} ({name})" if name else compound_id


def step_electron_fields(step: PlanStep) -> Dict[str, Any]:
    """把单步电子需求转换为路线汇总所需的标量字典。"""

    requirement = step.option.electron_requirement
    return {
        "electron_carrier_ids": ";".join(requirement.carrier_ids),
        "electron_requirement_classes": ";".join(requirement.requirement_classes),
        "electron_risk_level": requirement.risk_level,
        "electron_risk_score": requirement.risk_score,
        "electron_risk_evidence": "; ".join(requirement.evidence),
        "avoid_if_alternative_exists": requirement.avoid_if_alternative_exists,
        "requires_downstream_electron_design": requirement.requires_downstream_electron_design,
    }


def solution_burden_totals(solution: Solution) -> Dict[str, Any]:
    """汇总整条路线的辅因子、氧气和方向负担。"""

    return {
        "route_total_nadph_burden": sum(
            step.option.screening.nadph_burden for step in solution.steps
        ),
        "route_total_sam_burden": sum(
            step.option.screening.sam_burden for step in solution.steps
        ),
        "route_total_coa_burden": sum(
            step.option.screening.coa_burden for step in solution.steps
        ),
        "oxygen_required_steps": sum(
            1 for step in solution.steps if step.option.screening.oxygen_required
        ),
        "thermo_disfavored_steps": sum(
            1
            for step in solution.steps
            if step.option.screening.thermo_direction == "disfavored"
        ),
    }


def solution_electron_summary(solution: Solution) -> Dict[str, Any]:
    """汇总一条路线的电子传递系统风险。"""

    return summarize_solution_electron_requirements(
        step_electron_fields(step) for step in solution.steps
    )


def build_solution_summary_rows(
    result: SearchExecutionResult,
    target_compound: str,
    client: KeggRestClient,
) -> List[Dict[str, Any]]:
    """为每条候选路线生成一行总览，供 ``solutions.csv`` 使用。"""

    rows: List[Dict[str, Any]] = []

    for idx, solution in enumerate(result.solutions, start=1):
        heterologous_reactions = dedupe_keep_order(
            step.option.reaction.reaction_id
            for step in solution.steps
            if not step.is_endogenous
        )
        heterologous_kos = dedupe_keep_order(
            ko_id
            for step in solution.steps
            if not step.is_endogenous
            for ko_id in step.option.reaction.ko_ids
        )
        heterologous_ecs = dedupe_keep_order(
            ec
            for step in solution.steps
            if not step.is_endogenous
            for ec in step.option.reaction.enzyme_ecs
        )
        anchor_compounds = dedupe_keep_order(
            compound_id
            for step in solution.steps
            for compound_id in step.reachable_precursors
        )
        burden_totals = solution_burden_totals(solution)
        electron_summary = solution_electron_summary(solution)

        rows.append(
            {
                "solution_id": idx,
                "target_compound_id": target_compound,
                "target_compound_name": client.get_compound_name(target_compound),
                "total_steps": solution.total_steps,
                "heterologous_steps": solution.heterologous_steps,
                "heterologous_reaction_ids": ";".join(heterologous_reactions),
                "heterologous_ko_ids": ";".join(heterologous_kos),
                "heterologous_enzyme_ecs": ";".join(heterologous_ecs),
                "reaction_resolution_status": solution.reaction_resolution_status,
                "normalization_event_count": len(solution.normalization_events),
                "normalization_events": ";".join(solution.normalization_events),
                "blocking_reaction_count": len(solution.blocking_reaction_ids),
                "blocking_reaction_ids": ";".join(solution.blocking_reaction_ids),
                "eligible_for_recommendation": solution.reaction_ready,
                "reachable_anchor_compounds": ";".join(anchor_compounds),
                "reachable_anchor_labels": "; ".join(compound_label(kid, client) for kid in anchor_compounds),
                **burden_totals,
                "max_electron_risk_level": electron_summary["max_electron_risk_level"],
                "max_electron_risk_score": electron_summary["max_electron_risk_score"],
                "electron_system_status": electron_summary["electron_system_status"],
                "requires_downstream_electron_design": electron_summary[
                    "requires_downstream_electron_design"
                ],
            }
        )

    return rows


def build_solution_step_rows(
    solution_id: int,
    solution: Solution,
    client: KeggRestClient,
) -> List[Dict[str, Any]]:
    """将路线展开为 GEM 验证、酶选择和人工审阅所需的精简步骤行。"""

    rows: List[Dict[str, Any]] = []

    for step_index, step in enumerate(solution.steps, start=1):
        option = step.option
        reaction = option.reaction
        screening = option.screening
        row = {
            "solution_id": solution_id,
            "step_index": step_index,
            "status": "endogenous" if step.is_endogenous else "heterologous",
            "produced_compound_id": option.produced_compound,
            "produced_compound_name": client.get_compound_name(option.produced_compound),
            "reaction_id": reaction.reaction_id,
            "reaction_name": reaction.name,
            "equation": reaction.equation,
            "direction": option.direction,
            "oxygen_required": screening.oxygen_required,
            "thermo_direction": screening.thermo_direction,
            "screening_rule_hits": screening_hits_to_text(screening.rule_hits),
            "precursor_compound_ids": ";".join(option.precursor_compounds),
            "precursor_compound_labels": "; ".join(
                compound_label(kid, client) for kid in option.precursor_compounds
            ),
            "ko_ids": ";".join(reaction.ko_ids),
            "enzyme_ecs": ";".join(reaction.enzyme_ecs),
            "source_reaction_ids": ";".join(
                step.source_reaction_ids or (reaction.reaction_id,)
            ),
            "resolution_action": step.resolution_action,
            "resolution_evidence": ";".join(step.resolution_evidence),
        }
        # These source annotations are optional in KEGG.  Keep them when they
        # exist instead of manufacturing empty evidence in every exported row.
        if reaction.comment:
            row["reaction_comment"] = reaction.comment
        if reaction.module_ids:
            row["module_ids"] = ";".join(reaction.module_ids)
        if reaction.rhea_ids:
            row["rhea_ids"] = ";".join(reaction.rhea_ids)
        rows.append(row)

    return rows


def build_rejected_solution_rows(
    result: SearchExecutionResult,
) -> List[Dict[str, Any]]:
    """导出未通过反应规范化门禁的候选路线，供审计和修复使用。"""

    rows: List[Dict[str, Any]] = []
    seen: set[Tuple[str, ...]] = set()
    for solution in result.rejected_solutions:
        reaction_chain = tuple(
            step.option.reaction.reaction_id for step in solution.steps
        )
        if reaction_chain in seen:
            continue
        seen.add(reaction_chain)
        rows.append({
            "rejected_route_id": len(rows) + 1,
            "reaction_resolution_status": solution.reaction_resolution_status,
            "total_steps": solution.total_steps,
            "heterologous_steps": solution.heterologous_steps,
            "reaction_ids": ";".join(reaction_chain),
            "blocking_reaction_count": len(solution.blocking_reaction_ids),
            "blocking_reaction_ids": ";".join(solution.blocking_reaction_ids),
            "blocking_reasons": ";".join(solution.blocking_reasons),
            "normalization_events": ";".join(solution.normalization_events),
            "eligible_for_recommendation": False,
        })
    return rows


def build_electron_requirement_rows(result: SearchExecutionResult) -> List[Dict[str, Any]]:
    """生成仅包含高于零风险步骤的精简电子系统需求表。"""

    rows: List[Dict[str, Any]] = []
    for solution_id, solution in enumerate(result.solutions, start=1):
        for step_index, step in enumerate(solution.steps, start=1):
            requirement = step.option.electron_requirement
            if requirement.risk_score <= 0:
                continue
            rows.append(
                {
                    "solution_id": solution_id,
                    "step_index": step_index,
                    "reaction_id": step.option.reaction.reaction_id,
                    "electron_carrier_ids": ";".join(requirement.carrier_ids),
                    "electron_requirement_classes": ";".join(requirement.requirement_classes),
                    "electron_risk_level": requirement.risk_level,
                    "electron_risk_score": requirement.risk_score,
                    "electron_risk_evidence": "; ".join(requirement.evidence),
                    "requires_downstream_electron_design": (
                        requirement.requires_downstream_electron_design
                    ),
                }
            )
    return rows


def build_solution_electron_summary_rows(result: SearchExecutionResult) -> List[Dict[str, Any]]:
    """为每条路线生成独立于基础路线表的电子风险汇总。"""

    return [
        {
            "solution_id": solution_id,
            **solution_electron_summary(solution),
        }
        for solution_id, solution in enumerate(result.solutions, start=1)
    ]


def electron_search_stage(result: SearchExecutionResult) -> str:
    """由电子规避模式和兜底状态推导实际搜索阶段。"""

    if result.electron_avoidance_mode != "strict_with_fallback":
        return result.electron_avoidance_mode
    if result.electron_avoidance_fallback:
        return "prefer_expanded_fallback"
    return "strict_primary"


def write_outputs(
    target_compound: str,
    result: SearchExecutionResult,
    client: KeggRestClient,
    output_root: Path,
    run_args: Dict[str, Any],
) -> Path:
    """清理本目标的旧派生结果，并写出本轮搜索的 CSV/JSON 文件。

    清理范围仅限 ``kegg_gap_<target>`` 下已知的派生文件，以免旧路线、旧 GEM
    验证结果与本轮结果混用；KEGG 缓存和其他目标目录不受影响。
    """

    target_dir = safe_mkdir(output_root / f"kegg_gap_{target_compound}")
    # 同时清除旧版逐 solution 文件，避免它们与新的聚合表并存造成歧义。
    for pattern in (
        "solution_*_steps.csv",
        "solution_*_missing_enzymes.csv",
        "solutions.csv",
        "all_solution_steps.csv",
        "run_config.json",
        "report.md",
        "protein_candidates.csv",
        "route_electron_requirements.csv",
        "solution_electron_summary.csv",
        "rejected_reaction_routes.csv",
    ):
        for stale_path in target_dir.glob(pattern):
            if stale_path.is_file():
                stale_path.unlink()
    for validation_dir_name in ("gem_validation", "gem_validation_strict"):
        validation_dir = target_dir / validation_dir_name
        if not validation_dir.is_dir():
            continue
        for stale_path in validation_dir.glob("*.csv"):
            if stale_path.is_file():
                stale_path.unlink()

    # 总览与聚合步骤表分别服务路线选择、GEM 验证和后续酶筛选。
    summary_rows = build_solution_summary_rows(result, target_compound, client)
    pd.DataFrame(summary_rows).to_csv(target_dir / "solutions.csv", index=False, encoding="utf-8-sig")

    all_step_rows = [
        row
        for solution_id, solution in enumerate(result.solutions, start=1)
        for row in build_solution_step_rows(solution_id, solution, client)
    ]
    pd.DataFrame(all_step_rows).to_csv(target_dir / "all_solution_steps.csv", index=False, encoding="utf-8-sig")

    rejected_columns = [
        "rejected_route_id",
        "reaction_resolution_status",
        "total_steps",
        "heterologous_steps",
        "reaction_ids",
        "blocking_reaction_count",
        "blocking_reaction_ids",
        "blocking_reasons",
        "normalization_events",
        "eligible_for_recommendation",
    ]
    pd.DataFrame(
        build_rejected_solution_rows(result),
        columns=rejected_columns,
    ).to_csv(
        target_dir / "rejected_reaction_routes.csv",
        index=False,
        encoding="utf-8-sig",
    )

    # 单独导出存在电子风险的步骤，供后续电子传递伙伴设计直接消费。
    electron_requirement_rows = build_electron_requirement_rows(result)
    pd.DataFrame(electron_requirement_rows).to_csv(
        target_dir / "route_electron_requirements.csv",
        index=False,
        encoding="utf-8-sig",
    )
    electron_summary_rows = build_solution_electron_summary_rows(result)
    pd.DataFrame(electron_summary_rows).to_csv(
        target_dir / "solution_electron_summary.csv",
        index=False,
        encoding="utf-8-sig",
    )
    # 保存输入参数和实际兜底路径，保证搜索结果可复现、可审计。
    run_config = dict(run_args)
    run_config.update(
        {
            "search_mode_used": result.search_mode_used,
            "did_fallback": result.did_fallback,
            "electron_avoidance_mode": result.electron_avoidance_mode,
            "electron_avoidance_fallback": result.electron_avoidance_fallback,
            "electron_fallback_parameters": result.electron_fallback_parameters,
            "electron_search_stage_used": electron_search_stage(result),
            "screening_rule_version": SCREENING_RULE_VERSION,
            "reaction_resolution_version": REACTION_RESOLUTION_VERSION,
            "reaction_resolution_mode": result.reaction_resolution_mode,
            "rejected_reaction_route_count": len(
                build_rejected_solution_rows(result)
            ),
            "cycle_policy": CYCLE_POLICY,
            "target_compound_module_ids": list(result.target_compound_module_ids),
            "matched_target_module_ids": list(result.matched_target_module_ids),
            "expanded_module_ids": list(result.module_chain_context.expanded_module_ids),
            "module_chain_boundary_compounds": list(
                result.module_chain_context.boundary_compound_ids
            ),
            "allowed_module_reaction_count": len(
                result.module_chain_context.allowed_reaction_ids
            ),
            "direct_producing_reaction_ids": list(result.direct_producing_reaction_ids),
        }
    )
    (target_dir / "run_config.json").write_text(
        json.dumps(run_config, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    return target_dir


@tool
def kegg_gap_analyze(
    target: str,
    max_total_steps: int = 10,
    max_new_enzymes: int = MAX_HETEROLOGOUS_ENZYMES,
    max_solutions: int = 5,
    max_reactions_per_compound: int = 10,
    max_major_precursors: int = 4,
    max_routes_per_state: int = 3,
    max_module_chain_depth: int = 5,
    module_filter_mode: str = "prefer",
    electron_avoidance_mode: str = DEFAULT_ELECTRON_AVOIDANCE_MODE,
    reaction_resolution_mode: str = "strict",
):
    """
    分析目标 KEGG 化合物到宿主可产物的候选异源路径缺口。
    
    调用时机：用户给定目标化合物并要求通路/gap 分析。
    输入：target_compound_id、底盘/培养基/搜索参数。
    返回：ok、候选 solution 摘要、输出目录和后续蛋白选择所需文件。
    限制：长耗时计算；不写 manifest，用户确认方案后再调用 write_solution_to_manifest。
    """

    tool_name = "kegg_gap_analyze"
    monitor.report_start(tool_name, {"target": target})
    try:
        # 1. 参数校验
        if module_filter_mode not in {"prefer", "strict", "off"}:
            raise ValueError(
                f"module_filter_mode must be one of: prefer, strict, off. "
                f"Got: {module_filter_mode}"
            )
        electron_avoidance_mode = validate_electron_avoidance_mode(electron_avoidance_mode)
        reaction_resolution_mode = validate_reaction_resolution_mode(
            reaction_resolution_mode
        )

        target_compound = validate_target_compound_id(target)

        # 2. 统一转 Path
        model_path_obj = gem_model_file()
        output_root_obj = resolve_outputs_dir()
        reachable_path_obj = producible_kegg_compounds_file()
        cache_dir_obj = resolve_cache_dir()

        # 3. 加载输入
        monitor.report_running(tool_name, "正在加载可达化合物和内源 KEGG 反应...", progress=0.15)
        if not reachable_path_obj.exists() and (reachable_path_obj.parent / CHASSIS_ANALYSIS_RUNNING_MARKER).exists():
            monitor.report_running(
                tool_name,
                "正在等待底盘可达性分析完成...",
                progress=0.12,
            )
        wait_for_reachable_compounds_file(reachable_path_obj)
        reachable_compounds, _ = load_reachable_compounds(reachable_path_obj)
        endogenous_direction_index = load_endogenous_direction_index(
            model_path_obj,
            allowed_compartments=ENDOGENOUS_DIRECTION_COMPARTMENTS,
        )
        endogenous_reactions = set(endogenous_direction_index)

        # 4. 初始化 KEGG client
        client = KeggRestClient(cache_dir_obj)

        # 5. 搜索 KEGG gap solutions
        monitor.report_running(tool_name, "正在搜索候选异源合成路径...", progress=0.45)
        result = run_search_gap_analysis(
            target_compound=target_compound,
            reachable_compounds=reachable_compounds,
            endogenous_reactions=endogenous_reactions,
            client=client,
            max_total_steps=max_total_steps,
            max_new_enzymes=max_new_enzymes,
            max_solutions=max_solutions,
            max_reactions_per_compound=max_reactions_per_compound,
            max_major_precursors=max_major_precursors,
            ignored_common_compounds=set(IGNORED_COMMON_COMPOUNDS),
            module_filter_mode=module_filter_mode,
            max_routes_per_state=max_routes_per_state,
            max_module_chain_depth=max_module_chain_depth,
            electron_avoidance_mode=electron_avoidance_mode,
            reaction_resolution_mode=reaction_resolution_mode,
            endogenous_direction_index=endogenous_direction_index,
        )

        # 6. 手动构造 run_args，代替 vars(args)
        run_args: dict[str, Any] = {
            "target": target_compound,
            "model_path": str(model_path_obj),
            "reachable_file": str(reachable_path_obj),
            "max_total_steps": max_total_steps,
            "max_new_enzymes": max_new_enzymes,
            "max_solutions": max_solutions,
            "max_reactions_per_compound": max_reactions_per_compound,
            "max_major_precursors": max_major_precursors,
            "max_routes_per_state": max_routes_per_state,
            "max_module_chain_depth": max_module_chain_depth,
            "module_filter_mode": module_filter_mode,
            "electron_avoidance_mode": electron_avoidance_mode,
            "reaction_resolution_mode": reaction_resolution_mode,
            "endogenous_direction_mode": ENDOGENOUS_DIRECTION_MODE_GEM_BOUNDS,
            "endogenous_direction_compartments": sorted(
                ENDOGENOUS_DIRECTION_COMPARTMENTS
            ),
            "endogenous_direction_capable_reactions": direction_capable_reaction_count(
                endogenous_direction_index
            ),
        }

        # 7. 写输出文件
        monitor.report_running(tool_name, "正在写入 gap analysis 输出文件...", progress=0.9)
        output_dir = write_outputs(
            target_compound=target_compound,
            result=result,
            client=client,
            output_root=output_root_obj,
            run_args=run_args,
        )
        result_payload = {
            "target_compound": target_compound,
            "output_dir": str(output_dir.resolve()),
            "solutions_file": str((output_dir / "solutions.csv").resolve()),
            "all_solution_steps_file": str((output_dir / "all_solution_steps.csv").resolve()),
            "rejected_reaction_routes_file": str(
                (output_dir / "rejected_reaction_routes.csv").resolve()
            ),
            "solution_count": len(result.solutions),
            "search_mode_used": result.search_mode_used,
            "did_fallback": result.did_fallback,
            "electron_avoidance_mode": result.electron_avoidance_mode,
            "electron_avoidance_fallback": result.electron_avoidance_fallback,
            "reaction_resolution_mode": result.reaction_resolution_mode,
            "rejected_reaction_route_count": len(result.rejected_solutions),
        }
        monitor.report_end(tool_name, result_payload)
        return result_payload
    except Exception as exc:
        monitor.report_error(tool_name, exc)
        raise
