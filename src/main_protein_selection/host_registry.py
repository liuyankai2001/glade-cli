"""Validated, persistent host identities and read-only taxonomy snapshots."""
from __future__ import annotations

import hashlib
import json
import os
import tempfile
from pathlib import Path
from typing import Any

from filelock import FileLock, Timeout

REGISTRY_PATH = Path(__file__).resolve().parents[2] / 'data/chassis_taxonomy.json'
LOCK_DIRECTORY: Path | None = None
SCHEMA_VERSION = 'chassis_taxonomy.v1'


def _positive_id(value: Any) -> bool:
    return isinstance(value, int) and not isinstance(value, bool) and value > 0


def _name(value: Any) -> bool:
    return isinstance(value, str) and bool(value.strip())


def _unique_object(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f'duplicate registry JSON key: {key}')
        result[key] = value
    return result


def _validate_node(node):
    if (not isinstance(node, dict) or not _positive_id(node.get('taxonId'))
            or not _name(node.get('scientificName')) or not _name(node.get('rank'))):
        raise ValueError('invalid host registry taxonomy node')


def load_registry() -> dict:
    """Read fresh state; a missing or damaged file must never be silently reset."""
    path = Path(REGISTRY_PATH)
    try:
        payload = json.loads(path.read_text(encoding='utf-8-sig'), object_pairs_hook=_unique_object)
    except ValueError as exc:
        raise ValueError(f'invalid host registry {path}: {exc}') from exc
    if (not isinstance(payload, dict)
            or set(payload) != {'schema_version', 'chassis_taxon_presets', 'taxonomy_snapshots'}
            or payload['schema_version'] != SCHEMA_VERSION
            or not isinstance(payload['chassis_taxon_presets'], dict)
            or not isinstance(payload['taxonomy_snapshots'], dict)):
        raise ValueError(f'invalid host registry structure: {path}')
    ids = set()
    for key, record in payload['chassis_taxon_presets'].items():
        if (not _name(key) or key.strip() != key or not isinstance(record, dict)
                or set(record) != {'name', 'taxon_id'} or not _name(record['name'])
                or not _positive_id(record['taxon_id'])):
            raise ValueError(f'invalid host registry entry: {key}')
        if record['taxon_id'] in ids:
            raise ValueError(f'duplicate host registry taxon ID: {record["taxon_id"]}')
        ids.add(record['taxon_id'])
    for key, snapshot in payload['taxonomy_snapshots'].items():
        _validate_node(snapshot)
        if key != str(snapshot['taxonId']) or not isinstance(snapshot.get('lineage'), list):
            raise ValueError(f'invalid host registry taxonomy snapshot: {key}')
        for node in snapshot['lineage']:
            _validate_node(node)
    return payload


def get_host_preset(chassis_key: str) -> dict:
    preset = load_registry()['chassis_taxon_presets'].get(chassis_key)
    if preset is None:
        raise ValueError(f'Unknown chassis_key in host registry: {chassis_key}')
    return preset


def get_taxonomy_snapshot(taxon_id: int) -> dict | None:
    return load_registry()['taxonomy_snapshots'].get(str(taxon_id))


def _existing_key(payload: dict, taxon_id: int) -> str | None:
    return next((key for key, record in payload['chassis_taxon_presets'].items()
                 if record['taxon_id'] == taxon_id), None)


def _lock_path(path: Path) -> Path:
    identity = os.path.normcase(str(path.resolve())).encode('utf-8')
    directory = Path(LOCK_DIRECTORY) if LOCK_DIRECTORY is not None else Path(tempfile.gettempdir()) / 'glade-host-registry-locks'
    return directory / (hashlib.sha256(identity).hexdigest() + '.lock')


def _publish(path: Path, payload: dict) -> None:
    data = (json.dumps(payload, ensure_ascii=False, indent=2) + '\n').encode('utf-8')
    descriptor, temporary = tempfile.mkstemp(prefix='.chassis_taxonomy-', suffix='.tmp', dir=path.parent)
    try:
        with os.fdopen(descriptor, 'wb') as handle:
            handle.write(data)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def register_host(name: str, taxon_id: int) -> str:
    """Only candidate generation calls this mutating operation."""
    if not _name(name) or not _positive_id(taxon_id):
        raise ValueError('host name and positive integer taxon_id are required')
    path = Path(REGISTRY_PATH)
    payload = load_registry()
    existing = _existing_key(payload, taxon_id)
    if existing is not None:
        return existing
    lock = _lock_path(path)
    lock.parent.mkdir(parents=True, exist_ok=True)
    try:
        with FileLock(str(lock), timeout=30):
            payload = load_registry()
            existing = _existing_key(payload, taxon_id)
            if existing is not None:
                return existing
            key = f'taxon_{taxon_id}'
            if key in payload['chassis_taxon_presets']:
                raise ValueError(f'host registry key conflict: {key}')
            payload['chassis_taxon_presets'][key] = {'name': name.strip(), 'taxon_id': taxon_id}
            _publish(path, payload)
            return key
    except Timeout as exc:
        raise ValueError(f'host registry busy: {path}') from exc
