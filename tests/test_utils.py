from __future__ import annotations

import gzip

from PACKAGE.utils import smart_open


def test_smart_open_reads_real_gzip(tmp_path):
    path = tmp_path / "real.tsv.gz"
    with gzip.open(path, "wt") as f:
        f.write("read_id\tchrom\nfiber_a\tchr1\n")

    with smart_open(path) as f:
        assert f.readline() == "read_id\tchrom\n"
        assert f.readline() == "fiber_a\tchr1\n"


def test_smart_open_reads_plain_text_with_gz_suffix(tmp_path):
    path = tmp_path / "mislabeled.tsv.gz"
    path.write_text("read_id\tchrom\nfiber_a\tchr1\n")

    with smart_open(path) as f:
        assert f.readline() == "read_id\tchrom\n"
        assert f.readline() == "fiber_a\tchr1\n"
