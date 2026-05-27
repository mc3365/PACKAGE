"""Shared utilities used across PACKAGE modules.

Keep this module small. If something here grows beyond a few related helpers,
move it to its own module.
"""

from __future__ import annotations

import gzip
import logging
import sys
from pathlib import Path
from typing import IO


def get_logger(name: str = "PACKAGE") -> logging.Logger:
    """Get a configured logger.

    Idempotent: calling this multiple times with the same name returns the same logger
    without duplicating handlers.
    """
    logger = logging.getLogger(name)
    if logger.handlers:
        return logger
    handler = logging.StreamHandler(sys.stderr)
    handler.setFormatter(
        logging.Formatter("[%(asctime)s] %(name)s %(levelname)s: %(message)s", "%H:%M:%S")
    )
    logger.addHandler(handler)
    logger.setLevel(logging.INFO)
    return logger


def smart_open(filepath: str | Path, mode: str = "rt") -> IO:
    """Open a file with automatic gzip detection by extension.

    Args:
        filepath: Path to the file.
        mode: File open mode ('rt', 'rb', etc.).

    Returns:
        A file-like object. Caller is responsible for closing (use a with block).
    """
    filepath = Path(filepath)
    if filepath.suffix == ".gz":
        return gzip.open(filepath, mode)
    return open(filepath, mode)
