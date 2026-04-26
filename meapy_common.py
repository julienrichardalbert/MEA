#!/usr/bin/env python3
"""Common helpers for the MEApy transitional CLI."""

from __future__ import annotations

import shlex
import subprocess
from pathlib import Path
from typing import Iterable


REPO_ROOT = Path(__file__).resolve().parent


class MeapyError(RuntimeError):
    """Raised when a MEApy command cannot proceed."""


def require_path(path_text: str, kind: str = "file") -> Path:
    path = Path(path_text).expanduser().resolve()
    if kind == "file" and not path.is_file():
        raise MeapyError(f"Required file does not exist: {path}")
    if kind == "dir" and not path.is_dir():
        raise MeapyError(f"Required directory does not exist: {path}")
    return path


def run_command(cmd: Iterable[str]) -> None:
    cmd_list = list(cmd)
    if not cmd_list:
        raise MeapyError("Empty command provided.")
    print("[meapy] running:", " ".join(shlex.quote(x) for x in cmd_list))
    completed = subprocess.run(cmd_list, cwd=str(REPO_ROOT), check=False)
    if completed.returncode != 0:
        raise MeapyError(f"Command failed with code {completed.returncode}: {cmd_list[0]}")
