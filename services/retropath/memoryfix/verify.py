"""Check the generated plugin, agent and baked rule library using build provenance."""
import hashlib
import json
import os
from pathlib import Path


def sha(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main():
    base = Path(__file__).resolve().parent
    manifest = json.loads((base / "patch_manifest.json").read_text())
    paths = {"rdkit-nodes-memoryfix.jar": Path(manifest["original_plugin_path"]),
        "native-memory-agent.jar": base / "native-memory-agent.jar", "rr02.csv": Path("/opt/retropath/rules/rr02.csv")}
    expected = dict(manifest["artifacts_sha256"], **{"rr02.csv": manifest["rules_sha256"]})
    actual = {name: sha(path) for name, path in paths.items()}
    if actual != expected or os.environ.get("RETROPATH_RDKIT_PLUGIN_VERSION") != manifest["version"]:
        raise RuntimeError("Packaged memory repair, rules or runtime version differs from the build manifest")
    print(json.dumps({"verified": True, "version": manifest["version"], "sha256": actual}))


if __name__ == "__main__":
    main()
