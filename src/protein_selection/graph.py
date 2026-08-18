"""大肠杆菌酶蛋白补全工具的 LangGraph 主流程定义。"""

from typing import Literal, Protocol

from langgraph.graph import END, START, StateGraph
from langgraph.graph.state import CompiledStateGraph

from src.protein_selection.progress import emit_progress
from src.protein_selection.nodes.auxiliary_requirement import (
    FullUniProtLookup,
    StructuredOutputChatModel,
    build_auxiliary_requirement_node,
)
from src.protein_selection.nodes.main_research import (
    MainResearchAgentContextFactory,
    MainResearchAgentFactory,
    MainResearchAgentRunnable,
    build_main_research_node,
    return_reaction_mismatch,
)
from src.protein_selection.nodes.validation import (
    ReactionLookup,
    UniProtLookup,
    validate_input,
)
from src.protein_selection.state import ProteinSupplyState


# 节点名称集中定义，避免条件路由和节点注册使用不一致的字符串。
VALIDATION_NODE = "validate_input"
AUXILIARY_REQUIREMENT_NODE = "check_auxiliary_requirement"
RETURN_REACTION_MISMATCH_NODE = "return_reaction_mismatch"
MAIN_RESEARCH_NODE = "main_research_agent"
END_ROUTE = "end"


class UniProtWorkflowLookup(UniProtLookup, FullUniProtLookup, Protocol):
    """主图对 UniProt 客户端的完整能力要求。

    输入校验节点只需要精简记录，而辅助蛋白需求判断节点需要完整 UniProt
    注释，因此注入主图的客户端必须同时实现这两个查询协议。
    """


ValidationRoute = Literal[
    "check_auxiliary_requirement",
    "return_reaction_mismatch",
    "end",
]
RequirementRoute = Literal[
    "return_reaction_mismatch",
    "main_research_agent",
    "end",
]


def route_after_validation(state: ProteinSupplyState) -> ValidationRoute:
    """输入校验成功后继续分析，否则结束本次工作流。

    已由官方 Rhea 与 EC 映射确定为不匹配的输入直接进入确定性终态，不再分析
    辅助蛋白。其余有效输入继续进行亚基分析。

    ``invalid`` 和 ``service_error`` 都不会进入后续 LLM/智能体流程：前者表示
    输入本身无效，后者表示数据库服务异常。两种状态均由 CLI 根据图状态生成
    对应的错误信息和退出码。
    """

    if state.get("validation_status") == "valid":
        if state.get("preliminary_reaction_match") == "mismatched":
            emit_progress(
                "routing",
                "completed",
                "Rhea 与 EC 映射均明确冲突，进入反应不匹配快速路径",
            )
            return RETURN_REACTION_MISMATCH_NODE
        emit_progress("routing", "info", "输入有效，继续进行亚基分析")
        return AUXILIARY_REQUIREMENT_NODE
    emit_progress("routing", "skipped", "输入校验未通过，结束工作流")
    return END_ROUTE


def route_after_auxiliary_requirement(
    state: ProteinSupplyState,
) -> RequirementRoute:
    """根据初步辅助蛋白证据选择不匹配或完整研究路径。

    - ``mismatched``：若上游已确定 Rhea 与 EC 均冲突，直接输出不匹配终态。
    - ``required``：明确需要其他蛋白，交给研究智能体查找具体组分。
    - ``enhancing``：可独立催化但存在增效伙伴，交给研究智能体核验具体蛋白及
      其在大肠杆菌中的可用性。
    - ``undetermined``：现有注释不足，仍交给研究智能体补充数据库和文献证据。
    - ``not_required``：即使由旧状态或外部调用传入，也必须进入完整研究，不能
      仅凭亚基组成跳过辅助蛋白证据核验。
    - 缺失或未知状态：保守结束，避免把异常状态误判为成功。
    """

    if state.get("preliminary_reaction_match") == "mismatched":
        emit_progress(
            "routing",
            "completed",
            "已确认反应不匹配，跳过证据研究流程",
        )
        return RETURN_REACTION_MISMATCH_NODE

    status = state.get("requirement_status")
    if status in {"required", "enhancing", "undetermined", "not_required"}:
        emit_progress(
            "routing",
            "info",
            f"初步状态={status}，进入完整辅助蛋白证据研究流程",
        )
        return MAIN_RESEARCH_NODE
    emit_progress("routing", "error", "辅助需求状态无效，结束工作流")
    return END_ROUTE


def build_protein_supply_graph(
    *,
    llm: StructuredOutputChatModel,
    uniprot_client: UniProtWorkflowLookup,
    kegg_client: ReactionLookup,
    research_agent: MainResearchAgentRunnable | None = None,
    research_agent_factory: MainResearchAgentFactory | None = None,
    research_agent_context_factory: (
        MainResearchAgentContextFactory | None
    ) = None,
) -> CompiledStateGraph:
    """构建并编译输入校验、需求判断和证据研究主流程。"""

    def validation_node(state: ProteinSupplyState) -> ProteinSupplyState:
        """把带有外部数据库客户端依赖的校验函数适配为图节点。"""

        return validate_input(
            state,
            uniprot_client=uniprot_client,
            kegg_client=kegg_client,
        )

    # 使用统一的可序列化状态承载所有节点的输入和输出。
    graph = StateGraph(ProteinSupplyState)

    # 第一阶段：确定性校验 UniProt/KEGG 标识符及其数据库记录。
    graph.add_node(VALIDATION_NODE, validation_node)

    # 第二阶段：读取完整 UniProt 注释，仅判断是否需要另一种蛋白，不在这里
    # 搜索具体辅助蛋白。
    graph.add_node(
        AUXILIARY_REQUIREMENT_NODE,
        build_auxiliary_requirement_node(llm, uniprot_client),
    )

    # 确定性不匹配路径：Rhea 与 EC 均明确不相交时直接输出完整裁决。
    graph.add_node(
        RETURN_REACTION_MISMATCH_NODE,
        return_reaction_mismatch,
    )

    # 深度研究路径：由主管智能体协调数据库、文献、网页和宿主兼容性研究员，
    # 最终输出有证据支持的完整蛋白列表或未确定状态。
    graph.add_node(
        MAIN_RESEARCH_NODE,
        build_main_research_node(
            research_agent,
            research_agent_factory=research_agent_factory,
            research_agent_context_factory=research_agent_context_factory,
        ),
    )

    # 所有请求都从输入校验开始；校验失败时直接进入 END。
    graph.add_edge(START, VALIDATION_NODE)
    graph.add_conditional_edges(
        VALIDATION_NODE,
        route_after_validation,
        {
            AUXILIARY_REQUIREMENT_NODE: AUXILIARY_REQUIREMENT_NODE,
            RETURN_REACTION_MISMATCH_NODE: RETURN_REACTION_MISMATCH_NODE,
            END_ROUTE: END,
        },
    )

    # 需求判断完成后进行核心分流：除确定性反应不匹配外，全部进入完整研究；
    # 亚基组成不能单独证明不存在电子传递、成熟、装配或瞬时结合伙伴。
    graph.add_conditional_edges(
        AUXILIARY_REQUIREMENT_NODE,
        route_after_auxiliary_requirement,
        {
            RETURN_REACTION_MISMATCH_NODE: RETURN_REACTION_MISMATCH_NODE,
            MAIN_RESEARCH_NODE: MAIN_RESEARCH_NODE,
            END_ROUTE: END,
        },
    )

    # 两条业务分支都在写入最终状态后结束，不再合流到额外处理节点。
    graph.add_edge(RETURN_REACTION_MISMATCH_NODE, END)
    graph.add_edge(MAIN_RESEARCH_NODE, END)

    # compile() 会校验节点和边并生成可通过 invoke/ainvoke 执行的主图。
    return graph.compile()
