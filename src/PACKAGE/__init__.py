"""Compatibility namespace for the former PACKAGE project name.

New code should import :mod:`mei_fiber`. This namespace remains available so
existing analysis scripts can migrate without an immediate breaking change.
"""

from __future__ import annotations

import mei_fiber as _mei_fiber
from mei_fiber import FiberDatabase, __version__

# Let imports such as ``PACKAGE.db`` resolve to the renamed implementation.
__path__ = _mei_fiber.__path__

__all__ = ["FiberDatabase", "__version__"]
