"""Summarize the built distributions and check a tagged release against its index."""

import hashlib
import json
import os
from pathlib import Path


def check(files, published, limit):
    existing = {entry["filename"]: entry for entry in published}
    occupied = sum(entry["size"] for entry in existing.values())
    additional = 0
    for path in files:
        entry = existing.get(path.name)
        if entry is None:
            additional += path.stat().st_size
            continue
        with path.open("rb") as stream:
            digest = hashlib.file_digest(stream, "sha256").hexdigest()
        if digest != entry["hashes"]["sha256"]:
            raise ValueError(f"{path.name} already exists on the index with different contents")
    if additional and occupied + additional > limit:
        raise ValueError(f"The index holds {occupied:,} bytes and the release adds "
                         f"{additional:,}, exceeding the {limit:,}-byte project limit")
    return occupied, additional


if __name__ == "__main__":
    files = sorted(path for path in Path("dist").iterdir() if path.is_file())
    if not files:
        raise SystemExit("No built distributions found in dist")
    size = sum(path.stat().st_size for path in files)
    lines = [f"Built {len(files)} distributions totaling {size:,} bytes "
             f"({size / 2**20:,.1f} MiB)."]
    if os.environ.get("INDEX_HOST"):
        catalog = json.loads(Path(os.environ["INDEX_JSON"]).read_bytes())
        setting = ("TESTPYPI_PROJECT_LIMIT_BYTES" if os.environ["INDEX_HOST"] == "test.pypi.org"
                   else "PYPI_PROJECT_LIMIT_BYTES")
        limit = int(os.environ.get(setting) or 10 * 2**30)
        occupied, additional = check(files, catalog["files"], limit)
        lines.append(f"{os.environ['INDEX_HOST']} holds {occupied:,} bytes. "
                     f"This release adds {additional:,} bytes against a {limit:,}-byte limit.")
    summary = "\n\n".join(lines)
    print(summary)
    if os.environ.get("GITHUB_STEP_SUMMARY"):
        with open(os.environ["GITHUB_STEP_SUMMARY"], "a") as stream:
            print(summary, file=stream)
