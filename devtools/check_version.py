#!/usr/bin/env python3
"""Guard against version-string mistakes in the release workflow.

Two invariants are checked, selected by ``--mode``:

``dev`` (default; run on master and pull requests)
    ``master`` must always carry a development version, never a plain ``X.Y.Z``
    release version, once development is under way after a release. This catches
    the failure mode where the post-release ``.devN`` bump documented in the
    release procedure is skipped: master then reports an already-released (or
    otherwise final) version, so every source/development install misreports
    itself. This happened once, with master left on the released version for
    hundreds of commits.

    Concretely: if any commits exist after the most recent ``vX.Y.Z`` tag and
    the declared version is a plain ``X.Y.Z`` (no ``.devN``/pre-release suffix),
    the check fails. The transient release commit itself -- the tagged commit,
    with no commits after it -- is allowed to carry a plain version, so the
    check is independent of whether the branch or the tag is pushed first.

``release`` (run on tag pushes)
    The declared version must exactly match the tag being released. Because the
    version must equal the tag, a forgotten ``.dev0`` suffix (version
    ``0.9.14.dev0`` vs tag ``v0.9.14``) is rejected, while a deliberate
    pre-release tag (e.g. ``v0.9.7.b0``) whose version matches is allowed.

Exit status is 0 on success and 1 on any detected problem, so it can gate CI.
"""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
from pathlib import Path

VERSION_FILE = Path(__file__).resolve().parent.parent / "cnvlib" / "_version.py"

# ``__version__ = "X.Y.Z..."`` in cnvlib/_version.py.
VERSION_ASSIGN_RE = re.compile(r"""__version__\s*=\s*["']([^"']+)["']""")
# A plain final release: three dot-separated numbers, no dev/pre/post/local part.
RELEASE_VERSION_RE = re.compile(r"\d+\.\d+\.\d+")


def read_declared_version() -> str:
    """Return ``__version__`` from ``cnvlib/_version.py`` without importing it."""
    match = VERSION_ASSIGN_RE.search(VERSION_FILE.read_text())
    if not match:
        sys.exit(f"ERROR: no __version__ found in {VERSION_FILE}")
    return match.group(1)


def git(*args: str) -> str:
    """Run a git command and return its stripped stdout ('' on failure)."""
    try:
        out = subprocess.run(
            ["git", *args],
            capture_output=True,
            text=True,
            check=True,
        )
    except (subprocess.CalledProcessError, FileNotFoundError):
        return ""
    return out.stdout.strip()


def latest_release_tag() -> str | None:
    """Most recent ``vX.Y.Z`` tag reachable from HEAD, or None if none exist."""
    tag = git("describe", "--tags", "--abbrev=0", "--match", "v[0-9]*")
    return tag or None


def version_from_tag(tag: str) -> str:
    """Strip the leading ``v`` from a release tag (e.g. 'v0.9.14' -> '0.9.14')."""
    return tag.removeprefix("v")


def is_release_version(version: str) -> bool:
    """True if the version is a plain final release (``X.Y.Z``, no suffix)."""
    return RELEASE_VERSION_RE.fullmatch(version) is not None


def check_dev(version: str) -> int:
    tag = latest_release_tag()
    if tag is None:
        print(f"OK: no release tag yet; version is {version!r}.")
        return 0

    commits_after = git("rev-list", "--count", f"{tag}..HEAD")
    if commits_after in ("", "0"):
        # The tagged release commit itself (or git unavailable): a plain
        # release version is expected here, so nothing to enforce.
        print(f"OK: at release tag {tag!r} (no post-release commits).")
        return 0

    if is_release_version(version):
        if version == version_from_tag(tag):
            detail = f"which is already released as tag {tag!r}"
        else:
            detail = f"a plain release version (latest release tag is {tag!r})"
        print(
            f"ERROR: version is {version!r}, {detail}, but {commits_after} "
            "commit(s) have landed since the last release.\n"
            "       master must carry a development version -- bump to the next "
            "point release with a '.dev0' suffix (e.g. '0.9.15.dev0') per the "
            "release procedure.",
            file=sys.stderr,
        )
        return 1

    print(f"OK: version {version!r} is a development version.")
    return 0


def check_release(version: str, expected_tag: str) -> int:
    expected_version = version_from_tag(expected_tag)
    if version != expected_version:
        print(
            f"ERROR: declared version {version!r} does not match tag "
            f"{expected_tag!r} (expected {expected_version!r}).",
            file=sys.stderr,
        )
        return 1
    print(f"OK: release version {version!r} matches tag {expected_tag!r}.")
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--mode",
        choices=("dev", "release"),
        default="dev",
        help="Which invariant to check (default: dev).",
    )
    parser.add_argument(
        "--tag",
        help="For --mode=release: the tag being released (e.g. v0.9.14).",
    )
    args = parser.parse_args(argv)

    version = read_declared_version()
    if args.mode == "release":
        if not args.tag:
            parser.error("--mode=release requires --tag")
        return check_release(version, args.tag)
    return check_dev(version)


if __name__ == "__main__":
    raise SystemExit(main())
