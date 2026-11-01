"""Inspect supported geometry-backed basin registries without modifying them."""
import argparse
import json
from pyar.selection.basin_memory import _read_payload, BasinMemoryError


def main(argv=None, *, prog=None):
    parser = argparse.ArgumentParser(prog=prog or 'pyar basin-memory', description=__doc__)
    parser.add_argument('registry', nargs='+', help='supported basin_registry.json file(s)')
    args = parser.parse_args(argv)
    try:
        for path in args.registry:
            payload, _ = _read_payload(path)
            print(json.dumps(dict(path=path, schema_version=payload['schema_version'],
                                  entries=len(payload['entries']),
                                  geometry_entries=sum('geometry' in e for e in payload['entries'])), sort_keys=True))
    except BasinMemoryError as exc:
        parser.error(str(exc))
