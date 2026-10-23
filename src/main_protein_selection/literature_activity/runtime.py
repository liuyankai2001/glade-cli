"""Model configuration and literature tools for main-enzyme activity evidence."""

import os
from collections.abc import Mapping
from dataclasses import dataclass, replace
from pathlib import Path
from typing import cast

from dotenv import load_dotenv
from langchain.chat_models import init_chat_model
from langchain_core.language_models import BaseChatModel
from langchain_core.tools import BaseTool
from langchain_mcp_adapters.client import MultiServerMCPClient
from langchain_mcp_adapters.sessions import StdioConnection
from mcp.client.stdio import get_default_environment
from pydantic import SecretStr

PROJECT_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_ENV_PATH = PROJECT_ROOT / ".env"
SUPPORTED_MODEL_PROVIDER = "openai"
REQUIRED_ENV_KEYS = (
    "MODEL_PROVIDER",
    "AGENT_LLM_MODEL",
    "API_KEY",
    "BASE_URL",
)


class ConfigurationError(ValueError):
    """Raised when runtime configuration is missing or unsupported."""


@dataclass(frozen=True, slots=True)
class ModelSettings:
    """Validated model settings loaded from process environment or .env."""

    provider: str
    model: str
    api_key: SecretStr
    base_url: str


def load_model_settings(
    env_path: str | Path = DEFAULT_ENV_PATH,
) -> ModelSettings:
    """Load dotenv without overriding explicit process environment values."""

    load_dotenv(dotenv_path=env_path, override=False)
    values = {key: (os.getenv(key) or "").strip() for key in REQUIRED_ENV_KEYS}
    missing = [key for key, value in values.items() if not value]
    if missing:
        missing_text = ", ".join(missing)
        raise ConfigurationError(
            f"missing required model configuration: {missing_text}"
        )

    provider = values["MODEL_PROVIDER"].lower()
    if provider != SUPPORTED_MODEL_PROVIDER:
        raise ConfigurationError(
            "unsupported MODEL_PROVIDER: "
            f"{values['MODEL_PROVIDER']}; expected {SUPPORTED_MODEL_PROVIDER}"
        )

    return ModelSettings(
        provider=provider,
        model=values["AGENT_LLM_MODEL"],
        api_key=SecretStr(values["API_KEY"]),
        base_url=values["BASE_URL"].rstrip("/"),
    )


def build_chat_model(
    settings: ModelSettings,
    *,
    max_tokens: int = 4096,
    timeout_seconds: float = 90,
) -> BaseChatModel:
    """Build the configured provider through LangChain's model factory."""

    return init_chat_model(
        model=settings.model,
        model_provider=settings.provider,
        api_key=settings.api_key.get_secret_value(),
        base_url=settings.base_url,
        temperature=0,
        max_tokens=max_tokens,
        max_retries=0,
        timeout=timeout_seconds,
    )


TOOLUNIVERSE_SERVER_NAME = "bio-databases"
TOOLUNIVERSE_PACKAGE = "tooluniverse==1.4.1"
TOOLUNIVERSE_EXECUTABLE = "tooluniverse-smcp-stdio"

LITERATURE_RESEARCHER_TOOL_NAMES = (
    "PubMed_search_articles",
    "PubMed_get_article",
    "PubMed_get_related",
    "PubMed_get_cited_by",
    "EuropePMC_search_articles",
    "EuropePMC_get_fulltext_snippets",
    "EuropePMC_get_citations",
    "EuropePMC_get_references",
    "PubTator3_LiteratureSearch",
    "SemanticScholar_search_papers",
    "SemanticScholar_get_pdf_snippets",
    "Crossref_search_works",
    "Crossref_get_work",
    "Crossref_check_retraction",
)


@dataclass(frozen=True, slots=True)
class ToolUniverseConfig:
    """Runtime configuration for the isolated ToolUniverse MCP process."""

    uvx_command: str = "uvx"
    package: str = TOOLUNIVERSE_PACKAGE
    tool_names: tuple[str, ...] = LITERATURE_RESEARCHER_TOOL_NAMES
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


async def load_literature_researcher_tools(
    config: ToolUniverseConfig | None = None,
) -> list[BaseTool]:
    """Load and verify the complete literature tool allowlist."""

    literature_config = replace(
        config or ToolUniverseConfig(), tool_names=LITERATURE_RESEARCHER_TOOL_NAMES
    )
    available_tools = await load_tooluniverse_tools(literature_config)
    tools_by_name = {tool.name: tool for tool in available_tools}
    missing = [
        name for name in LITERATURE_RESEARCHER_TOOL_NAMES if name not in tools_by_name
    ]
    if missing:
        missing_names = ", ".join(missing)
        raise RuntimeError(
            f"ToolUniverse did not expose required literature tools: {missing_names}"
        )

    return [tools_by_name[name] for name in LITERATURE_RESEARCHER_TOOL_NAMES]
