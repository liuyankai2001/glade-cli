"""Load the local open-websearch MCP server as LangChain tools."""

from dataclasses import dataclass
import sys
from typing import Literal, cast

from langchain_core.tools import BaseTool
from langchain_mcp_adapters.client import MultiServerMCPClient
from langchain_mcp_adapters.sessions import StdioConnection


OPEN_WEBSEARCH_SERVER_NAME = "web-search"
OPEN_WEBSEARCH_NPM_PACKAGE = "open-websearch@2.1.11"


@dataclass(frozen=True, slots=True)
class OpenWebSearchConfig:
    """Runtime configuration for the local open-websearch MCP process."""

    npm_package: str = OPEN_WEBSEARCH_NPM_PACKAGE
    default_search_engine: str = "bing"
    allowed_search_engines: tuple[str, ...] | None = None
    proxy_url: str | None = None
    search_mode: Literal["request", "auto", "playwright"] = "auto"
    playwright_navigation_timeout_ms: int = 20000


def build_open_websearch_connection(
    config: OpenWebSearchConfig | None = None,
    *,
    platform: str | None = None,
) -> StdioConnection:
    """Build a cross-platform stdio connection for open-websearch.

    Windows uses ``cmd /c`` because ``npx`` is normally installed as
    ``npx.cmd``. ``-y`` prevents an interactive install prompt from blocking
    the MCP handshake.
    """

    settings = config or OpenWebSearchConfig()
    current_platform = platform or sys.platform
    env = {
        "MODE": "stdio",
        "DEFAULT_SEARCH_ENGINE": settings.default_search_engine,
        "SEARCH_MODE": settings.search_mode,
        "PLAYWRIGHT_NAVIGATION_TIMEOUT_MS": str(
            settings.playwright_navigation_timeout_ms
        ),
    }

    if settings.allowed_search_engines:
        env["ALLOWED_SEARCH_ENGINES"] = ",".join(
            settings.allowed_search_engines
        )
    if settings.proxy_url:
        env.update({"USE_PROXY": "true", "PROXY_URL": settings.proxy_url})

    if current_platform.startswith("win"):
        command = "cmd"
        args = ["/c", "npx", "-y", settings.npm_package]
    else:
        command = "npx"
        args = ["-y", settings.npm_package]

    return cast(
        StdioConnection,
        {
            "transport": "stdio",
            "command": command,
            "args": args,
            "env": env,
            "encoding": "utf-8",
        },
    )


def build_open_websearch_client(
    config: OpenWebSearchConfig | None = None,
) -> MultiServerMCPClient:
    """Create the LangChain MCP client without starting the subprocess."""

    return MultiServerMCPClient(
        {
            OPEN_WEBSEARCH_SERVER_NAME: build_open_websearch_connection(config),
        },
        handle_tool_errors=True,
    )


async def load_open_websearch_tools(
    config: OpenWebSearchConfig | None = None,
) -> list[BaseTool]:
    """Start open-websearch, discover its tools, and adapt them for LangChain.

    The adapter opens a fresh stdio MCP session when an adapted tool is called,
    so the returned tools remain usable after this discovery call finishes.
    """

    client = build_open_websearch_client(config)
    return await client.get_tools(server_name=OPEN_WEBSEARCH_SERVER_NAME)
