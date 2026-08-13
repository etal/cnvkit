#!/usr/bin/env python3
"""Check that every ``cnvkit.py`` command in the manual is one argparse accepts.

The manual's runnable examples are the part of the documentation readers copy,
and the part that silently rots: an option can be renamed, gain an argument or
move to another subcommand without anything in ``doc/`` noticing. This script
extracts every ``cnvkit.py ...`` invocation from ``doc/*.rst`` and feeds it to
the real parser from :mod:`cnvlib.commands`, so the manual is validated against
the program rather than against a copy of its option list. It found
``--smooth-bootstrap`` documented without its argument six months after the
option stopped being a flag.

The check is deliberately shallow. It asks only whether argparse accepts the
command line, not whether the command does what the surrounding prose claims;
files named in the examples do not exist and are never opened. What it defends
is the mechanical layer -- option spellings, arities, ``choices`` values and
subcommand names -- which is exactly the layer prose cannot keep in step by
hand.

Exit status is 0 when every invocation parses and 1 otherwise, so it can gate
CI. It runs in the documentation environment, ahead of ``sphinx-build``, because
that is the one CI job triggered by changes to both ``doc/`` and ``cnvlib/``.
"""

from __future__ import annotations

import argparse
import contextlib
import io
import itertools
import re
import shlex
import sys
from pathlib import Path

from cnvlib.commands import AP

#: A line that starts a command, optionally behind a shell prompt.
INVOCATION_RE = re.compile(r"^\s*\$?\s*(cnvkit\.py\s.*)$")

#: Shell operators that end the CNVkit command and begin something else.
TAIL_RE = re.compile(r"\s(?:\||>|>>|&&|;)\s|\s(?:\||>|>>)")

#: A brace expansion such as ``Sample.cn{s,r}``, which the examples use to name
#: a .cns and .cnr pair in one word.
BRACE_RE = re.compile(r"\{([^{}]*,[^{}]*)\}")


def iter_doc_lines(doc_dir: Path):
    """Yield ``(path, line number, text)`` with backslash continuations joined."""
    for path in sorted(doc_dir.glob("*.rst")):
        lines = path.read_text().splitlines()
        index = 0
        while index < len(lines):
            text, last = lines[index], index
            while text.rstrip().endswith("\\") and last + 1 < len(lines):
                last += 1
                text = text.rstrip()[:-1] + " " + lines[last].strip()
            yield path, index + 1, text
            index = last + 1


def expand_braces(word: str) -> list[str]:
    """Expand one level of ``{a,b}`` alternation, as the shell would.

    ``Sample.cn{s,r}`` becomes two words, which matters: the examples rely on
    the expansion to supply both the ``-s`` segment file and the positional
    ``.cnr`` file.
    """
    match = BRACE_RE.search(word)
    if not match:
        return [word]
    prefix, suffix = word[: match.start()], word[match.end() :]
    return list(
        itertools.chain.from_iterable(
            expand_braces(prefix + alternative + suffix)
            for alternative in match.group(1).split(",")
        )
    )


def to_argv(command: str) -> list[str] | None:
    """Turn a documented command line into an argv list, or None to skip it.

    Three conventions of the manual are normalized away. A command may be
    piped into another program, in which case only the CNVkit part is ours to
    check; a fragment may be bracketed to mark it optional, as in the
    ``batch`` pipeline listing; and file pairs are written with brace
    expansion. Templated invocations, which contain placeholders rather than
    arguments, are skipped: the WDL examples in ``doc/docker.rst`` are not
    shell commands.
    """
    if "${" in command:
        return None
    tail = TAIL_RE.search(command)
    if tail:
        command = command[: tail.start()]
    try:
        words = shlex.split(command)
    except ValueError:
        return None
    words = [word.strip("[]") for word in words if word.strip("[]")]
    return list(itertools.chain.from_iterable(expand_braces(w) for w in words))


def parse_error(parser: argparse.ArgumentParser, argv: list[str]) -> str | None:
    """Return argparse's complaint about `argv`, or None if it is accepted.

    ``--help`` exits zero before parsing the rest, which the manual does on
    purpose when it shows readers how to get usage text; that is a pass.

    argparse builds every message around ``prog``, which it takes from
    ``sys.argv[0]`` -- this checker -- and bakes into each subparser at import
    time. Reporting an error against ``check_doc_cli.py segmetrics`` would name
    a command the reader cannot run, so the prefix is rewritten to the real one.
    """
    stderr = io.StringIO()
    try:
        with (
            contextlib.redirect_stderr(stderr),
            contextlib.redirect_stdout(io.StringIO()),
        ):
            parser.parse_args(argv[1:])
    except SystemExit as exit_:
        if exit_.code:
            message = stderr.getvalue().strip().splitlines()
            if not message:
                return f"exit status {exit_.code}"
            return message[-1].replace(parser.prog, "cnvkit.py", 1)
    return None


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument(
        "doc_dir",
        nargs="?",
        default=Path(__file__).resolve().parent.parent / "doc",
        type=Path,
        help="Directory of .rst files to check. [Default: the doc/ directory]",
    )
    args = ap.parse_args(argv)

    checked = 0
    failures = []
    for path, lineno, text in iter_doc_lines(args.doc_dir):
        match = INVOCATION_RE.match(text)
        if not match:
            continue
        command = match.group(1)
        parsed = to_argv(command)
        if parsed is None:
            continue
        checked += 1
        error = parse_error(AP, parsed)
        if error:
            failures.append((path, lineno, command, error))

    for path, lineno, command, error in failures:
        print(f"{path}:{lineno}\n    {command}\n    {error}", file=sys.stderr)
    print(f"checked {checked} documented invocations; {len(failures)} rejected")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
