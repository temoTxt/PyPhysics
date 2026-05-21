"""Validation tools: validate_wikilinks, regenerate_indexes, regenerate_tracker."""

import subprocess

from pyphysics_mcp.repo import repo_root


def _run_script(script_rel_path: str, args: list[str] | None = None) -> dict:
    """Shell out to a repo Python script via `uv run python`."""
    args = args or []
    script = repo_root() / script_rel_path
    if not script.exists():
        return {"error": f"script not found: {script_rel_path}"}
    cmd = ["uv", "run", "python", str(script), *args]
    res = subprocess.run(cmd, capture_output=True, text=True, cwd=repo_root())
    return {
        "status": "ok" if res.returncode == 0 else "error",
        "returncode": res.returncode,
        "stdout": res.stdout.strip(),
        "stderr": res.stderr.strip(),
    }


def validate_wikilinks(paths: list[str] | None = None) -> dict:
    """Run the wikilink validator. Pass --scan for each requested path.

    Returns: {status, returncode, stdout, stderr, broken_count}.
    """
    args: list[str] = []
    for p in paths or []:
        args.extend(["--scan", p])
    out = _run_script("Roadmapping/History/_tools/validate_wikilinks.py", args)
    # Extract broken-link count from stderr if present.
    broken = 0
    for line in out.get("stderr", "").splitlines():
        if "broken wikilink" in line:
            try:
                broken = int(line.split()[0])
            except (ValueError, IndexError):
                pass
            break
    out["broken_count"] = broken
    return out


def regenerate_indexes() -> dict:
    """Run build_dataview_indexes.py to refresh the three index pages."""
    return _run_script("Roadmapping/History/_tools/build_dataview_indexes.py")


def regenerate_tracker() -> dict:
    """Run update_acquisition_tracker.py to refresh Historical_Papers/Acquisition_Tracker.md."""
    return _run_script("Roadmapping/History/Bibliography/update_acquisition_tracker.py")
