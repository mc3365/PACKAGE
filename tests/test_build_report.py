"""Tests for full-build benchmark reporting."""

from __future__ import annotations

import json

from PACKAGE.benchmark.build_report import parse_gnu_time, write_build_report


def test_parse_gnu_time(tmp_path):
    resource_log = tmp_path / "time.txt"
    resource_log.write_text(
        "\tUser time (seconds): 12.34\n"
        "\tMaximum resident set size (kbytes): 987654\n"
        "\tElapsed (wall clock) time (h:mm:ss or m:ss): 1:02.50\n"
        "\tExit status: 0\n"
    )

    result = parse_gnu_time(resource_log)

    assert result["user_time_(seconds)"] == 12.34
    assert result["maximum_resident_set_size_(kbytes)"] == 987654
    assert result["elapsed_(wall_clock)_time_(h:mm:ss_or_m:ss)"] == "1:02.50"
    assert result["exit_status"] == 0


def test_write_build_report(tiny_db_path, tmp_path):
    config = tmp_path / "benchmark.yaml"
    config.write_text(
        f"""output_dir: {tiny_db_path.parent}
output_file: {tiny_db_path.name}
samples:
  - name: d0
    layers:
      nucleosomes: {tmp_path / 'unused.csv'}
"""
    )
    resource_log = tmp_path / "time.txt"
    resource_log.write_text("\tMaximum resident set size (kbytes): 12345\n")
    output = tmp_path / "report.json"

    report = write_build_report(
        config,
        resource_log,
        output,
        started_at="2026-06-22T10:00:00-04:00",
        finished_at="2026-06-22T10:01:00-04:00",
        elapsed_seconds=60,
        git_commit="abc123",
    )

    assert report["samples"] == {"d0": 2}
    assert report["gnu_time"]["maximum_resident_set_size_(kbytes)"] == 12345
    assert json.loads(output.read_text())["git_commit"] == "abc123"
