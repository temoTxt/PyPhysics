"""Animation tools: render_scene."""

import subprocess
from pathlib import Path

from pyphysics_mcp.repo import animations_dir, relpath


def render_scene(scene_name: str, quality: str = "ql") -> dict:
    """Render a Manim scene from Roadmapping/Animations/manim_scenes/.

    Args:
        scene_name: scene filename stem OR ClassName. Examples:
            'hist_faraday_induction' (filename)
            'HistFaradayInduction' (class name; we'll try to locate the file)
        quality: 'ql' (480p, fast), 'qh' (1080p, slow), 'qm' (720p).

    Returns: {status, scene_file, class_name, command, returncode, video_path?, stdout, stderr}.
    """
    if quality not in ("ql", "qm", "qh"):
        return {"error": f"quality must be 'ql', 'qm', or 'qh'; got {quality!r}"}

    scenes_dir = animations_dir() / "manim_scenes"
    if not scenes_dir.is_dir():
        return {"error": f"manim_scenes directory not found at {scenes_dir}"}

    # Resolve scene file + class name.
    scene_file: Path | None = None
    class_name: str | None = None

    # Try as filename stem first.
    candidate = scenes_dir / f"{scene_name}.py"
    if candidate.exists():
        scene_file = candidate
        class_name = _extract_class_name(candidate)
    else:
        # Try as class name across all scene files.
        for sf in scenes_dir.glob("*.py"):
            extracted = _extract_class_name(sf)
            if extracted == scene_name:
                scene_file = sf
                class_name = extracted
                break

    if scene_file is None or class_name is None:
        return {"error": f"scene {scene_name!r} not found under {relpath(scenes_dir)}"}

    cmd = [
        "uv", "run", "manim", f"-{quality}",
        "--media_dir", "rendered",
        str(scene_file.relative_to(animations_dir())),
        class_name,
    ]
    res = subprocess.run(cmd, capture_output=True, text=True, cwd=animations_dir())

    # Locate the rendered video.
    quality_dir = {"ql": "480p15", "qm": "720p30", "qh": "1080p60"}[quality]
    video = animations_dir() / "rendered" / "videos" / scene_file.stem / quality_dir / f"{class_name}.mp4"

    return {
        "status": "ok" if res.returncode == 0 else "error",
        "scene_file": relpath(scene_file),
        "class_name": class_name,
        "command": " ".join(cmd),
        "returncode": res.returncode,
        "video_path": relpath(video) if video.exists() else None,
        "stdout": res.stdout.strip()[-2000:],  # tail-truncate
        "stderr": res.stderr.strip()[-2000:],
    }


def _extract_class_name(scene_file: Path) -> str | None:
    """Extract the first `class Foo(Scene):` name from a scene file."""
    import re
    text = scene_file.read_text(encoding="utf-8")
    m = re.search(r"^class\s+(\w+)\s*\(.*Scene\s*\)\s*:", text, re.MULTILINE)
    return m.group(1) if m else None
