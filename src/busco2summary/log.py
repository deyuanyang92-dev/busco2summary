"""Minimal logging to stderr plus an optional run log file."""

from __future__ import annotations

import sys
import time
from pathlib import Path
from typing import Optional, TextIO

_LOG: Optional[TextIO] = None
QUIET = False


def set_logfile(path: Optional[Path]) -> None:
    global _LOG
    if _LOG:
        _LOG.close()
    _LOG = open(path, "a", encoding="utf-8") if path else None


def _emit(level: str, msg: str) -> None:
    line = f"[{time.strftime('%H:%M:%S')}] {level:<5} {msg}"
    if not QUIET or level in ("WARN", "ERROR"):
        print(line, file=sys.stderr, flush=True)
    if _LOG:
        _LOG.write(line + "\n")
        _LOG.flush()


def info(msg: str) -> None:
    _emit("INFO", msg)


def warn(msg: str) -> None:
    _emit("WARN", msg)


def error(msg: str) -> None:
    _emit("ERROR", msg)
