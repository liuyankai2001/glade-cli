"""Environment-backed runtime configuration for Protein Supply."""

from dataclasses import dataclass
import os
from pathlib import Path

from dotenv import load_dotenv
from langchain.chat_models import init_chat_model
from langchain_core.language_models import BaseChatModel
from pydantic import SecretStr


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_ENV_PATH = PROJECT_ROOT / ".env"
DEFAULT_CACHE_DIR = PROJECT_ROOT / ".cache" / "protein_supply"
CACHE_DIR_ENV_KEY = "CACHE_DIR"
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


def load_cache_dir(
    env_path: str | Path = DEFAULT_ENV_PATH,
) -> Path:
    """Load the cache root from ``CACHE_DIR`` in the process or dotenv.

    Relative values are resolved against the directory containing ``.env``
    so the configured location is stable even when the CLI is started from a
    different working directory. Explicit process environment values retain
    precedence over dotenv values.
    """

    dotenv_path = Path(env_path)
    load_dotenv(dotenv_path=dotenv_path, override=False)
    configured = (os.getenv(CACHE_DIR_ENV_KEY) or "").strip()
    if not configured:
        return DEFAULT_CACHE_DIR

    cache_dir = Path(os.path.expandvars(configured)).expanduser()
    if not cache_dir.is_absolute():
        cache_dir = dotenv_path.resolve().parent / cache_dir
    return cache_dir.resolve(strict=False)


def load_model_settings(
    env_path: str | Path = DEFAULT_ENV_PATH,
) -> ModelSettings:
    """Load dotenv without overriding explicit process environment values."""

    load_dotenv(dotenv_path=env_path, override=False)
    values = {
        key: (os.getenv(key) or "").strip()
        for key in REQUIRED_ENV_KEYS
    }
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
