"""Repo-path resolution + shared helpers (YAML frontmatter parsing).

Path resolution order:
  1. PYPHYSICS_REPO_PATH env var if set and contains Roadmapping/History/.
  2. Walk up from cwd looking for Roadmapping/History/Bibliography/.
  3. Error with a clear remediation message.
"""

import os
import re
from functools import lru_cache
from pathlib import Path

import yaml

FRONTMATTER_RE = re.compile(r"^---\n(.*?)\n---", re.DOTALL)


class RepoNotFound(RuntimeError):
    pass


@lru_cache(maxsize=1)
def repo_root() -> Path:
    env = os.environ.get("PYPHYSICS_REPO_PATH")
    if env:
        p = Path(env).expanduser().resolve()
        if (p / "Roadmapping" / "History" / "Bibliography").is_dir():
            return p
        raise RepoNotFound(
            f"PYPHYSICS_REPO_PATH={env!r} does not contain Roadmapping/History/Bibliography/"
        )

    cur = Path.cwd().resolve()
    for ancestor in (cur, *cur.parents):
        if (ancestor / "Roadmapping" / "History" / "Bibliography").is_dir():
            return ancestor

    raise RepoNotFound(
        "Could not locate the PyPhysics repo. Set PYPHYSICS_REPO_PATH=/path/to/repo "
        "or run pyphysics-mcp from inside the repo tree."
    )


def parse_frontmatter(md_path: Path) -> dict:
    """Return YAML frontmatter as a dict; {} if missing or unparseable."""
    text = md_path.read_text(encoding="utf-8")
    m = FRONTMATTER_RE.match(text)
    if not m:
        return {}
    try:
        return yaml.safe_load(m.group(1)) or {}
    except yaml.YAMLError:
        return {}


def split_frontmatter(md_path: Path) -> tuple[dict, str]:
    """Return (frontmatter_dict, body_str)."""
    text = md_path.read_text(encoding="utf-8")
    m = FRONTMATTER_RE.match(text)
    if not m:
        return {}, text
    try:
        meta = yaml.safe_load(m.group(1)) or {}
    except yaml.YAMLError:
        meta = {}
    body = text[m.end():].lstrip("\n")
    return meta, body


# Convenience path accessors.

def bib_dir() -> Path:
    return repo_root() / "Roadmapping" / "History" / "Bibliography"


def history_dir() -> Path:
    return repo_root() / "Roadmapping" / "History"


def forward_dir() -> Path:
    return repo_root() / "Roadmapping" / "History" / "Forward"


def animations_dir() -> Path:
    return repo_root() / "Roadmapping" / "Animations"


def relpath(p: Path) -> str:
    """Path relative to repo root (POSIX-style for JSON output)."""
    try:
        return p.resolve().relative_to(repo_root()).as_posix()
    except ValueError:
        return p.as_posix()
