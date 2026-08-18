"""Query-scoped persistent MCP sessions for the research workflow."""

from collections.abc import Sequence
from contextlib import AsyncExitStack
from dataclasses import dataclass, replace
from types import TracebackType

from langchain_core.tools import BaseTool
from langchain_mcp_adapters.tools import load_mcp_tools

from src.protein_selection.integrations.open_websearch import (
    OPEN_WEBSEARCH_SERVER_NAME,
    OpenWebSearchConfig,
    build_open_websearch_client,
)
from src.protein_selection.integrations.tooluniverse import (
    TOOLUNIVERSE_SERVER_NAME,
    ToolUniverseConfig,
    build_tooluniverse_client,
)


@dataclass(frozen=True, slots=True)
class ResearchToolNames:
    """Tool allowlists assigned to the four isolated researchers."""

    bio_database: tuple[str, ...]
    literature: tuple[str, ...]
    web: tuple[str, ...]
    host_compatibility: tuple[str, ...]

    def __post_init__(self) -> None:
        for group_name, names in (
            ("bio_database", self.bio_database),
            ("literature", self.literature),
            ("web", self.web),
            ("host_compatibility", self.host_compatibility),
        ):
            if not names:
                raise ValueError(f"{group_name} tool allowlist cannot be empty")
            if len(names) != len(set(names)):
                raise ValueError(
                    f"{group_name} tool allowlist contains duplicate names"
                )

    @property
    def tooluniverse_union(self) -> tuple[str, ...]:
        """Return a stable de-duplicated union for one server process."""

        return tuple(
            dict.fromkeys(
                (
                    *self.bio_database,
                    *self.literature,
                    *self.host_compatibility,
                )
            )
        )


@dataclass(frozen=True, slots=True)
class ResearchMCPTools:
    """Session-bound tools partitioned for each research subagent."""

    bio_database: tuple[BaseTool, ...]
    literature: tuple[BaseTool, ...]
    web: tuple[BaseTool, ...]
    host_compatibility: tuple[BaseTool, ...]


class ResearchMCPRuntime:
    """Own one ToolUniverse and one web MCP session for a single query."""

    def __init__(
        self,
        *,
        tool_names: ResearchToolNames,
        tooluniverse_config: ToolUniverseConfig | None = None,
        web_config: OpenWebSearchConfig | None = None,
    ) -> None:
        self._tool_names = tool_names
        self._tooluniverse_config = tooluniverse_config
        self._web_config = web_config
        self._exit_stack: AsyncExitStack | None = None
        self._tools: ResearchMCPTools | None = None

    @property
    def tools(self) -> ResearchMCPTools:
        """Return active tools and reject access outside the session scope."""

        if self._tools is None:
            raise RuntimeError("research MCP runtime is not active")
        return self._tools

    async def __aenter__(self) -> ResearchMCPTools:
        if self._exit_stack is not None:
            raise RuntimeError("research MCP runtime is already active")

        exit_stack = AsyncExitStack()
        try:
            base_tooluniverse_config = (
                self._tooluniverse_config or ToolUniverseConfig()
            )
            shared_tooluniverse_config = replace(
                base_tooluniverse_config,
                tool_names=self._tool_names.tooluniverse_union,
            )
            tooluniverse_client = build_tooluniverse_client(
                shared_tooluniverse_config
            )
            tooluniverse_session = await exit_stack.enter_async_context(
                tooluniverse_client.session(TOOLUNIVERSE_SERVER_NAME)
            )
            tooluniverse_tools = await load_mcp_tools(
                tooluniverse_session,
                server_name=TOOLUNIVERSE_SERVER_NAME,
                handle_tool_errors=True,
            )

            web_client = build_open_websearch_client(self._web_config)
            web_session = await exit_stack.enter_async_context(
                web_client.session(OPEN_WEBSEARCH_SERVER_NAME)
            )
            web_tools = await load_mcp_tools(
                web_session,
                server_name=OPEN_WEBSEARCH_SERVER_NAME,
                handle_tool_errors=True,
            )

            tools = ResearchMCPTools(
                bio_database=_select_tools(
                    tooluniverse_tools,
                    self._tool_names.bio_database,
                    group_name="bio database",
                ),
                literature=_select_tools(
                    tooluniverse_tools,
                    self._tool_names.literature,
                    group_name="literature",
                ),
                web=_select_tools(
                    web_tools,
                    self._tool_names.web,
                    group_name="web",
                ),
                host_compatibility=_select_tools(
                    tooluniverse_tools,
                    self._tool_names.host_compatibility,
                    group_name="host compatibility",
                ),
            )
        except BaseException:
            await exit_stack.aclose()
            raise

        self._exit_stack = exit_stack
        self._tools = tools
        return tools

    async def __aexit__(
        self,
        exc_type: type[BaseException] | None,
        exc: BaseException | None,
        traceback: TracebackType | None,
    ) -> bool:
        exit_stack = self._exit_stack
        self._exit_stack = None
        self._tools = None
        if exit_stack is None:
            return False
        return await exit_stack.__aexit__(exc_type, exc, traceback)


def _select_tools(
    available_tools: Sequence[BaseTool],
    required_names: Sequence[str],
    *,
    group_name: str,
) -> tuple[BaseTool, ...]:
    tools_by_name = {tool.name: tool for tool in available_tools}
    missing = [name for name in required_names if name not in tools_by_name]
    if missing:
        missing_names = ", ".join(missing)
        raise RuntimeError(
            f"MCP server did not expose required {group_name} tools: "
            f"{missing_names}"
        )
    return tuple(tools_by_name[name] for name in required_names)
