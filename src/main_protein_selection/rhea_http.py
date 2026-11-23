"""Direct Rhea transport, kept separate from proxy-aware protein/KEGG sessions."""
from __future__ import annotations

import time
import requests

from .settings import GLADE_CONTACT_EMAIL_ENV, RHEA_HTTP_CONFIG, env_value


def new_rhea_session() -> requests.Session:
    session = requests.Session()
    session.trust_env = False
    contact = env_value(GLADE_CONTACT_EMAIL_ENV)
    session.headers["User-Agent"] = f"GLADE/2.3 ({contact})" if contact else "GLADE/2.3"
    return session


def request_response(session, url, *, params=None, timeout=30, retries=3, stop_on_client_error=True):
    last_error = None
    for attempt in range(max(1, retries)):
        try:
            response = session.get(url, params=params, timeout=timeout)
            response.raise_for_status()
            if response.status_code != 200:
                raise requests.HTTPError(f"Unexpected HTTP status {response.status_code}", response=response)
            return response
        except requests.RequestException as exc:
            last_error = exc
            status = exc.response.status_code if exc.response is not None else None
            if stop_on_client_error and status is not None and 400 <= status < 500 and status != 429:
                break
            if attempt + 1 < retries:
                time.sleep(RHEA_HTTP_CONFIG.sleep_seconds * (2 ** attempt))
    raise requests.RequestException(f"Rhea/official source request failed: {last_error}") from last_error


def reject_html(text: str) -> None:
    lowered = text.lstrip("\ufeff \r\n\t").lower()
    if lowered.startswith("<") or "<html" in lowered or "cf_chl_opt" in lowered:
        raise ValueError("HTML/browser challenge returned instead of scientific data")
