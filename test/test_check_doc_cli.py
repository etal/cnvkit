#!/usr/bin/env python
"""Tests for devtools/check_doc_cli.py, the documented-command guard.

The script exists because a documented command line can stop working without
anything in ``doc/`` changing: ``--smooth-bootstrap`` became an integer option
in 9fc3456 and ``doc/pipeline.rst`` went on showing it as a bare flag for six
months. The tests below defend two halves of that contract.

  * The normalizer accepts the conventions the manual actually uses -- shell
    pipelines, bracketed optional fragments, brace-expanded file pairs, prompts
    -- because a checker that rejects valid prose would be turned off.
  * The check rejects the drift classes it exists to catch, and specifically
    the arity change that motivated it. A guard that cannot fail is worse than
    no guard, since it reports success either way.

``check_doc_cli`` lives in ``devtools/``, which is not an importable package,
so it is loaded by file path as ``test_check_version`` does.
"""

import importlib.util
from pathlib import Path

import pytest

_SCRIPT = Path(__file__).resolve().parent.parent / "devtools" / "check_doc_cli.py"
_spec = importlib.util.spec_from_file_location("check_doc_cli", _SCRIPT)
cdc = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(cdc)


def write_doc(tmp_path, command):
    """Write a one-command .rst file and return its directory."""
    (tmp_path / "example.rst").write_text(f"Example\n=======\n\n::\n\n    {command}\n")
    return tmp_path


@pytest.mark.parametrize(
    ("command", "expected"),
    [
        # A pipeline: only the CNVkit part is ours to check.
        (
            "cnvkit.py breaks S.cnr S.cns | cut -f1 > genes.txt",
            ["cnvkit.py", "breaks", "S.cnr", "S.cns"],
        ),
        # Brace expansion supplies both the -s segment file and the positional
        # .cnr file from one word, so it must expand rather than be stripped.
        (
            "cnvkit.py scatter -s Sample.cn{s,r}",
            ["cnvkit.py", "scatter", "-s", "Sample.cns", "Sample.cnr"],
        ),
        # Brackets mark an optional fragment in the batch pipeline listing.
        (
            "cnvkit.py autobin *.bam [--annotate refFlat.txt]",
            ["cnvkit.py", "autobin", "*.bam", "--annotate", "refFlat.txt"],
        ),
        # A shell prompt is stripped by the line pattern, not by to_argv.
        ("cnvkit.py version", ["cnvkit.py", "version"]),
    ],
)
def test_to_argv_normalizes_manual_conventions(command, expected):
    assert cdc.to_argv(command) == expected


def test_to_argv_skips_templated_invocations():
    """WDL placeholders are not shell words and cannot be parsed as arguments."""
    assert cdc.to_argv("cnvkit.py batch ${bam_files} -r ${reference}") is None


def test_expand_braces_is_recursive():
    assert cdc.expand_braces("a{1,2}b{x,y}") == ["a1bx", "a1by", "a2bx", "a2by"]


def test_expand_braces_leaves_plain_words_alone():
    assert cdc.expand_braces("Sample.cnr") == ["Sample.cnr"]


def test_continuation_lines_are_joined(tmp_path):
    """A backslash-continued example is one command, not two fragments."""
    (tmp_path / "c.rst").write_text(
        "T\n=\n\n::\n\n    cnvkit.py batch S.bam \\\n        -r ref.cnn -d out\n"
    )
    texts = [text.strip() for _, _, text in cdc.iter_doc_lines(tmp_path)]
    assert "cnvkit.py batch S.bam  -r ref.cnn -d out" in texts


def test_accepts_the_current_manual():
    """The shipped documentation must pass; this is the check's real subject."""
    doc_dir = Path(__file__).resolve().parent.parent / "doc"
    assert cdc.main([str(doc_dir)]) == 0


def test_rejects_the_arity_change_that_motivated_the_check(tmp_path, capsys):
    """--smooth-bootstrap took no argument until 9fc3456; now it requires one."""
    doc_dir = write_doc(
        tmp_path, "cnvkit.py segmetrics S.cnr -s S.cns --ci --smooth-bootstrap -o o.cns"
    )
    assert cdc.main([str(doc_dir)]) == 1
    err = capsys.readouterr().err
    assert "--smooth-bootstrap" in err
    # The reader is told about cnvkit.py, not about the checker that ran.
    assert "cnvkit.py segmetrics: error:" in err


@pytest.mark.parametrize(
    "command",
    [
        "cnvkit.py batch *.bam --normals Normal.bam",  # renamed option
        "cnvkit.py segment S.cnr --nonexistent",  # removed option
        "cnvkit.py fix a.cnn b.cnn -r r.cnn --alpha 0.05",  # option on wrong command
        "cnvkit.py call S.cns --method quantile",  # value outside choices
        "cnvkit.py genemetric S.cnr",  # misspelled subcommand
        "cnvkit.py export vcf4 S.cns",  # misspelled nested subcommand
    ],
)
def test_rejects_drift(tmp_path, command):
    assert cdc.main([str(write_doc(tmp_path, command))]) == 1


@pytest.mark.parametrize(
    "command",
    [
        "cnvkit.py call S.cns --method threshold -t=-1,0,1",
        "cnvkit.py export vcf S.cns -o S.vcf",
        "cnvkit.py scatter -s Sample.cn{s,r} -c chr7 -g BRAF",
        "cnvkit.py batch -h",  # the manual shows readers how to get usage text
    ],
)
def test_accepts_valid_commands(tmp_path, command):
    assert cdc.main([str(write_doc(tmp_path, command))]) == 0
