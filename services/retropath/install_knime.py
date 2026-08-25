from __future__ import annotations

import argparse
import hashlib
import os
import subprocess
import tarfile
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import requests


KNIME_URL = (
    "https://zenodo.org/api/records/7564938/files/"
    "knime_4.7.0.linux.gtk.x86_64.tar.gz/content"
)
KNIME_FILENAME = "knime_4.7.0.linux.gtk.x86_64.tar.gz"
KNIME_SIZE = 590_209_872
KNIME_MD5 = "9f2bee3e470182a1003e79dfc1653090"
KNIME_TOP_LEVEL_DIR = "knime_4.7.0"
P2_REPOSITORIES = (
    "https://update.knime.com/analytics-platform/4.7/",
    "https://update.knime.com/community-contributions/trusted/4.7/",
)
P2_FEATURES = (
    "org.knime.features.chem.types.feature.group/4.7.1.v202301311311",
    "org.knime.features.datageneration.feature.group/4.7.0.v202211082353",
    "org.knime.features.python.feature.group/4.7.1.v202301311311",
    "org.rdkit.knime.feature.feature.group/4.9.1.v202312081930",
)


def md5_file(path: Path) -> str:
    digest = hashlib.md5()  # noqa: S324 - verifies the upstream Zenodo record
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _download_range(
    start: int,
    end: int,
    destination: Path,
    *,
    attempts: int = 20,
) -> None:
    expected_size = end - start + 1
    for attempt in range(1, attempts + 1):
        current_size = destination.stat().st_size if destination.exists() else 0
        if current_size == expected_size:
            return
        if current_size > expected_size:
            destination.unlink()
            current_size = 0
        range_start = start + current_size
        try:
            with requests.get(
                KNIME_URL,
                headers={"Range": f"bytes={range_start}-{end}"},
                stream=True,
                allow_redirects=True,
                timeout=(30, 120),
            ) as response:
                if response.status_code != 206:
                    raise RuntimeError(
                        f"Zenodo did not honor byte range {range_start}-{end}: "
                        f"HTTP {response.status_code}"
                    )
                with destination.open("ab") as output:
                    for chunk in response.iter_content(chunk_size=1024 * 1024):
                        if chunk:
                            output.write(chunk)
            current_size = destination.stat().st_size
            if current_size == expected_size:
                print(
                    f"KNIME range complete: {start}-{end} ({expected_size} bytes)",
                    flush=True,
                )
                return
        except (requests.RequestException, OSError, RuntimeError) as exc:
            current_size = destination.stat().st_size if destination.exists() else 0
            print(
                f"KNIME range {start}-{end} attempt {attempt}/{attempts} "
                f"interrupted at {current_size}/{expected_size} bytes: {exc}",
                flush=True,
            )
        if attempt < attempts:
            time.sleep(min(5 * attempt, 30))
    current_size = destination.stat().st_size if destination.exists() else 0
    raise RuntimeError(
        f"KNIME range {start}-{end} incomplete: {current_size}/{expected_size} bytes"
    )


def download_with_resume(destination: Path, *, workers: int = 8) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    prefix_size = destination.stat().st_size if destination.exists() else 0
    if prefix_size == KNIME_SIZE:
        return
    if prefix_size > KNIME_SIZE:
        destination.unlink()
        prefix_size = 0

    remaining = KNIME_SIZE - prefix_size
    chunk_size = (remaining + workers - 1) // workers
    ranges: list[tuple[int, int, Path]] = []
    for index in range(workers):
        start = prefix_size + index * chunk_size
        if start >= KNIME_SIZE:
            break
        end = min(start + chunk_size - 1, KNIME_SIZE - 1)
        part = destination.parent / f"{KNIME_FILENAME}.range-{start}-{end}"
        ranges.append((start, end, part))

    print(
        f"KNIME parallel download: existing prefix={prefix_size}, "
        f"ranges={len(ranges)}",
        flush=True,
    )
    with ThreadPoolExecutor(max_workers=len(ranges)) as executor:
        futures = {
            executor.submit(_download_range, start, end, part): (start, end)
            for start, end, part in ranges
        }
        for future in as_completed(futures):
            future.result()

    assembled = destination.parent / f"{KNIME_FILENAME}.assembled"
    assembled.unlink(missing_ok=True)
    with assembled.open("wb") as output:
        if prefix_size:
            with destination.open("rb") as prefix:
                for chunk in iter(lambda: prefix.read(1024 * 1024), b""):
                    output.write(chunk)
        for _, _, part in ranges:
            with part.open("rb") as source:
                for chunk in iter(lambda: source.read(1024 * 1024), b""):
                    output.write(chunk)
    os.replace(assembled, destination)
    for _, _, part in ranges:
        part.unlink(missing_ok=True)
    size = destination.stat().st_size
    if size != KNIME_SIZE:
        raise RuntimeError(
            f"assembled KNIME archive has unexpected size: {size}/{KNIME_SIZE} bytes"
        )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cache-dir", type=Path, required=True)
    parser.add_argument("--install-dir", type=Path, required=True)
    args = parser.parse_args()

    archive = args.cache_dir / KNIME_FILENAME
    download_with_resume(archive)
    actual_md5 = md5_file(archive)
    if actual_md5 != KNIME_MD5:
        archive.unlink(missing_ok=True)
        raise RuntimeError(
            f"KNIME MD5 mismatch: expected {KNIME_MD5}, got {actual_md5}"
        )
    print(f"KNIME archive verified: md5={actual_md5}", flush=True)

    args.install_dir.mkdir(parents=True, exist_ok=True)
    knime_root = args.install_dir / KNIME_TOP_LEVEL_DIR
    if knime_root.exists():
        raise RuntimeError(f"KNIME install directory is not empty: {knime_root}")
    with tarfile.open(archive, "r:gz") as package:
        top_level = {
            Path(member.name).parts[0]
            for member in package.getmembers()
            if member.name and Path(member.name).parts
        }
        if top_level != {KNIME_TOP_LEVEL_DIR}:
            raise RuntimeError(
                f"unexpected KNIME archive top-level entries: {sorted(top_level)}"
            )
        package.extractall(args.install_dir)

    knime_executable = knime_root / "knime"
    p2_dir = knime_root / "p2"
    if not knime_executable.is_file() or not p2_dir.is_dir():
        raise RuntimeError("extracted KNIME package is missing its executable or p2 directory")

    command = [
        str(knime_executable),
        "-nosplash",
        "-consoleLog",
        "-application",
        "org.eclipse.equinox.p2.director",
        "-repository",
        ",".join(P2_REPOSITORIES),
        "-bundlepool",
        str(p2_dir),
        "-destination",
        str(knime_root),
        "-i",
        ",".join(P2_FEATURES),
    ]
    print("Installing pinned KNIME and RDKit features", flush=True)
    subprocess.run(command, check=True)

    required = {
        "RDKit nodes": list(
            (knime_root / "p2" / "plugins").glob("org.rdkit.knime.nodes_*.jar")
        ),
        "KNIME chemistry feature": list(
            (knime_root / "p2" / "features").glob(
                "org.knime.features.chem.types_*"
            )
        ),
    }
    missing = [name for name, matches in required.items() if not matches]
    if missing:
        raise RuntimeError(f"KNIME plugin verification failed: {', '.join(missing)}")
    print("KNIME plugin verification succeeded", flush=True)


if __name__ == "__main__":
    main()
