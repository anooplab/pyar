"""Bind CI distribution evidence to a candidate; never authorizes publishing."""

import argparse
import hashlib
import json
from pathlib import Path


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def validate(record, sha, version, directory, required):
    if (
        record["schema_version"] != 1
        or record["sha"] != sha
        or record["version"] != version
    ):
        raise ValueError("Candidate identity mismatch")
    if not required or set(record["checks"]) != set(required):
        raise ValueError("Required check inventory mismatch")
    for check in record["checks"].values():
        if check != {"sha": sha, "conclusion": "success"}:
            raise ValueError("Failed, stale or incomplete required check")
    paths = sorted(Path(directory).glob("*.whl")) + sorted(
        Path(directory).glob("*.tar.gz")
    )
    if len(paths) != 2 or {p.suffix for p in paths} != {".whl", ".gz"}:
        raise ValueError("Exactly one wheel and one sdist required")
    if record["artifacts"] != {p.name: digest(p) for p in paths}:
        raise ValueError("Artifact inventory/hash mismatch")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("--sha", required=True)
    parser.add_argument("--version", required=True)
    parser.add_argument("--results", type=Path, required=True, help="GitHub needs JSON")
    parser.add_argument("--required", nargs="+", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    results = json.loads(args.results.read_text())
    record = {
        "schema_version": 1,
        "sha": args.sha,
        "version": args.version,
        "checks": {
            name: {"sha": args.sha, "conclusion": job["result"]}
            for name, job in results.items()
        },
        "artifacts": {
            p.name: digest(p) for p in args.directory.iterdir() if p.is_file()
        },
    }
    validate(record, args.sha, args.version, args.directory, args.required)
    args.output.write_text(json.dumps(record, indent=2) + "\n")


if __name__ == "__main__":
    main()
