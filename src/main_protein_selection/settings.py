"""Service defaults used by the standalone main-enzyme selector."""

from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path

from dotenv import dotenv_values


@dataclass(frozen=True, slots=True)
class HttpConfig:
    timeout_seconds: float = 30.0
    retries: int = 3
    sleep_seconds: float = 0.2


UNIPROT_HTTP_CONFIG = HttpConfig()
RHEA_HTTP_CONFIG = HttpConfig()
KEGG_HTTP_CONFIG = HttpConfig()
SELENZYME_HTTP_CONFIG = HttpConfig(timeout_seconds=60.0)

UNIPROT_SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search"
UNIPROT_PAGE_SIZE = 500
RHEA_REST_URL = "https://www.rhea-db.org/rhea"
KEGG_REST_BASE_URL = "https://rest.kegg.jp"

GLADE_CONTACT_EMAIL_ENV = "GLADE_CONTACT_EMAIL"
SELENZYME_REST_URL_ENV = "SELENZYME_REST_URL"
DEFAULT_ENV_PATH = Path(__file__).resolve().parents[2] / ".env"


def env_value(name: str) -> str:
    """Return an environment value, falling back to the repository ``.env``."""

    process_value = str(os.getenv(name) or "").strip()
    if process_value:
        return process_value
    if not DEFAULT_ENV_PATH.is_file():
        return ""
    return str(dotenv_values(DEFAULT_ENV_PATH).get(name) or "").strip()


def get_selenzyme_rest_url() -> str:
    """Return the configured Selenzyme endpoint or fail at fallback time."""

    value = env_value(SELENZYME_REST_URL_ENV)
    if not value:
        raise RuntimeError(
            f"{SELENZYME_REST_URL_ENV} is required for Selenzyme fallback"
        )
    return value


__all__ = [
    "DEFAULT_ENV_PATH",
    "GLADE_CONTACT_EMAIL_ENV",
    "HttpConfig",
    "KEGG_HTTP_CONFIG",
    "KEGG_REST_BASE_URL",
    "RHEA_HTTP_CONFIG",
    "RHEA_REST_URL",
    "SELENZYME_HTTP_CONFIG",
    "SELENZYME_REST_URL_ENV",
    "UNIPROT_HTTP_CONFIG",
    "UNIPROT_PAGE_SIZE",
    "UNIPROT_SEARCH_URL",
    "env_value",
    "get_selenzyme_rest_url",
]
