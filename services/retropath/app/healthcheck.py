from __future__ import annotations

import json
import sys
import urllib.request


def main() -> int:
    try:
        with urllib.request.urlopen("http://127.0.0.1:8000/health", timeout=5) as response:
            payload = json.load(response)
    except Exception:
        return 1
    return 0 if response.status == 200 and payload.get("ready") is True else 1


if __name__ == "__main__":
    sys.exit(main())

