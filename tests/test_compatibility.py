"""Compatibility checks for the former PACKAGE namespace."""

from __future__ import annotations


def test_legacy_package_namespace_reexports_current_api():
    from mei_fiber import FiberDatabase, __version__
    from PACKAGE import FiberDatabase as LegacyFiberDatabase
    from PACKAGE import __version__ as legacy_version

    assert LegacyFiberDatabase is FiberDatabase
    assert legacy_version == __version__ == "0.5.0"


def test_legacy_nested_import_resolves_to_current_implementation():
    from PACKAGE.db import FiberDatabase as LegacyFiberDatabase

    from mei_fiber.db import FiberDatabase

    assert LegacyFiberDatabase is FiberDatabase
