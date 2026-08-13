"""Render every ``cnvkit.py`` subcommand's ``--help`` output into the manual.

The ``cli-reference`` directive walks the argparse tree built in
:mod:`cnvlib.commands` and emits one section per command, each containing that
command's help text exactly as the terminal shows it.

Help text is inserted as a literal block, i.e. as an opaque text node rather
than as reStructuredText. CNVkit's help strings were written for an 80-column
terminal, not for publication: they contain glob patterns (``*.cnr``), option
names (``--diagram``) and quoted values with apostrophes, all of which docutils
would otherwise interpret as markup. Parsing them produces broken emphasis, en
dashes in place of option prefixes, and spurious definition lists -- silently,
in the case of the third-party extensions surveyed before this module was
written. A literal block cannot be misread, so the page is correct by
construction rather than by escaping every help string.
"""

from __future__ import annotations

import argparse
import os
from typing import TYPE_CHECKING

from docutils import nodes
from docutils.parsers.rst import Directive

from cnvlib.commands import AP

if TYPE_CHECKING:
    from collections.abc import Iterator


def format_help(parser: argparse.ArgumentParser, prog: str) -> str:
    """Return `parser`'s help text as it would appear under the name `prog`.

    Two details of argparse have to be overridden for the result to be
    reproducible. It derives ``prog`` from ``sys.argv[0]``, which during a
    documentation build is the build tool, and it wraps help text to the
    caller's terminal, which would otherwise make the rendered page depend on
    the environment that happened to build it. Both are restored afterward,
    since the parser is a module-level singleton shared with whatever else
    imports :mod:`cnvlib.commands`.
    """
    old_prog = parser.prog
    old_columns = os.environ.get("COLUMNS")
    parser.prog = prog
    os.environ["COLUMNS"] = "80"
    try:
        return parser.format_help()
    finally:
        parser.prog = old_prog
        if old_columns is None:
            del os.environ["COLUMNS"]
        else:
            os.environ["COLUMNS"] = old_columns


def walk_parsers(
    parser: argparse.ArgumentParser, prog: str
) -> Iterator[tuple[str, argparse.ArgumentParser]]:
    """Yield ``(program name, parser)`` for `parser` and every subcommand.

    Deprecated command aliases are registered by assigning an existing parser
    to a second key in the subparser map, so ``genemetrics`` and ``gainloss``
    are one object under two names. Yielding both would publish the deprecated
    spelling as though it were current, and -- because the two share a ``prog``
    attribute -- would leave whichever name came last in the rendered usage
    line of both. Identity, not name, decides what has already been seen.
    """
    yield prog, parser
    for action in parser._actions:
        if isinstance(action, argparse._SubParsersAction):
            seen: set[int] = set()
            for name, subparser in action.choices.items():
                if id(subparser) in seen:
                    continue
                seen.add(id(subparser))
                yield from walk_parsers(subparser, f"{prog} {name}")


class CliReference(Directive):
    """Emit a section of ``--help`` output for each CNVkit subcommand."""

    has_content = False

    def run(self) -> list[nodes.Node]:
        sections = []
        for prog, parser in walk_parsers(AP, "cnvkit.py"):
            help_text = format_help(parser, prog)
            section = nodes.section(ids=[nodes.make_id(prog)], names=[prog])
            section += nodes.title(text=prog)
            block = nodes.literal_block(help_text, help_text)
            # Without an explicit lexer Pygments guesses, and colors fragments
            # of the help text as though they were Python source.
            block["language"] = "none"
            section += block
            sections.append(section)
        return sections


def setup(app):
    app.add_directive("cli-reference", CliReference)
    return {"parallel_read_safe": True, "parallel_write_safe": True}
