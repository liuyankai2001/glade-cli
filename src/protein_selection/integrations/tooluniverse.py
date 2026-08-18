"""Load a restricted ToolUniverse MCP server as LangChain tools."""

from collections.abc import Mapping
from dataclasses import dataclass
from typing import cast

from langchain_core.tools import BaseTool
from langchain_mcp_adapters.client import MultiServerMCPClient
from langchain_mcp_adapters.sessions import StdioConnection
from mcp.client.stdio import get_default_environment


TOOLUNIVERSE_SERVER_NAME = "bio-databases"
TOOLUNIVERSE_PACKAGE = "tooluniverse==1.4.1"
TOOLUNIVERSE_EXECUTABLE = "tooluniverse-smcp-stdio"

# This explicit allowlist keeps the database researcher away from ToolUniverse's
# web, literature, code-execution, and agent-composition tools.
TOOLUNIVERSE_DEFAULT_TOOL_NAMES = (
    "UniProt_get_entry_by_accession",
    "UniProt_get_function_by_accession",
    "UniProt_get_organism_by_accession",
    "UniProt_get_sequence_by_accession",
    "UniProt_search",
    "proteins_api_search",
    "KEGG_link_entries",
    "KEGG_convert_ids",
    "KEGG_get_reaction",
    "KEGG_get_enzyme",
    "Rhea_get_reaction",
    "Rhea_get_reaction_participants",
    "Rhea_search_by_ec",
    "ComplexPortal_search_complexes",
    "ComplexPortal_get_complex",
    "intact_search_interactions",
    "intact_get_interactions",
    "STRING_map_identifiers",
    "STRING_get_protein_interactions",
    "InterPro_get_entries_for_protein",
    "InterPro_get_protein_domains",
)


@dataclass(frozen=True, slots=True)
class ToolUniverseConfig:
    """Runtime configuration for the isolated ToolUniverse MCP process."""

    uvx_command: str = "uvx"
    package: str = TOOLUNIVERSE_PACKAGE
    tool_names: tuple[str, ...] = TOOLUNIVERSE_DEFAULT_TOOL_NAMES
    extra_env: Mapping[str, str] | None = None


def build_tooluniverse_connection(
    config: ToolUniverseConfig | None = None,
) -> StdioConnection:
    """Build the stdio connection without starting ToolUniverse.

    ToolUniverse remains outside the project's Python environment and is
    launched through ``uvx`` at the pinned version. The stdio server disables
    hooks by default; ``--no-hooks`` is deliberately not passed because it is
    not a supported argument in ToolUniverse 1.4.1.
    """

    settings = config or ToolUniverseConfig()
    if not settings.tool_names:
        raise ValueError("ToolUniverse tool allowlist cannot be empty")

    env = get_default_environment()
    if settings.extra_env:
        env.update(settings.extra_env)
    env.update(
        {
            "PYTHONIOENCODING": "utf-8",
            "PYTHONUTF8": "1",
        }
    )

    args = [
        "--from",
        settings.package,
        TOOLUNIVERSE_EXECUTABLE,
        "--include-tools",
        *settings.tool_names,
        "--no-search",
    ]
    return cast(
        StdioConnection,
        {
            "transport": "stdio",
            "command": settings.uvx_command,
            "args": args,
            "env": env,
            "encoding": "utf-8",
            "encoding_error_handler": "strict",
        },
    )


def build_tooluniverse_client(
    config: ToolUniverseConfig | None = None,
) -> MultiServerMCPClient:
    """Create the LangChain MCP client without starting the subprocess."""

    return MultiServerMCPClient(
        {
            TOOLUNIVERSE_SERVER_NAME: build_tooluniverse_connection(config),
        },
        handle_tool_errors=True,
    )


async def load_tooluniverse_tools(
    config: ToolUniverseConfig | None = None,
) -> list[BaseTool]:
    """Start ToolUniverse and adapt its allowlisted tools for LangChain."""

    client = build_tooluniverse_client(config)
    return await client.get_tools(server_name=TOOLUNIVERSE_SERVER_NAME)
