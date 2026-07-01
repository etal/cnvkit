#!/usr/bin/env python
"""Tests for devtools/check_version.py, the release version-string guard.

This script gates CI against two real failure modes, one of which already bit
the project: master was left on the released version ``0.9.13`` for 276 commits
because the post-release ``.dev0`` bump was skipped, so every source install
misreported its version. The tests below defend that contract:

  * ``--mode dev`` detects a *stuck release* -- the declared version is a plain
    ``X.Y.Z`` release with commits landed after the latest tag -- and rejects
    it, whether the version equals that tag or has simply run ahead of it.
  * ``--mode dev`` allows the transient *release commit* (version == tag with
    nothing after it) and any development version, so the check is push-order
    independent.
  * ``--mode release`` treats the tag as the source of truth: it accepts any
    declared version equal to ``version_from_tag(tag)`` -- including a
    prerelease such as ``0.9.7.b0`` when the tag itself carries that suffix --
    and rejects anything that does not match.

``git(...)`` and ``read_declared_version()`` are monkeypatched, so no real
repository or subprocess is touched: the tests are hermetic and fast.
"""

import importlib.util
from pathlib import Path

import pytest

# check_version lives in devtools/, which is not an importable package, so load
# it by file path. The module resolves its own paths from __file__ at import
# time, so the autouse chdir fixture in conftest.py does not affect it.
_SCRIPT = Path(__file__).resolve().parent.parent / "devtools" / "check_version.py"
_spec = importlib.util.spec_from_file_location("check_version", _SCRIPT)
cv = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(cv)


def fake_git(describe="", rev_list=""):
    """Build a stand-in for cv.git that answers describe/rev-list by subcommand.

    ``describe`` is returned for ``git describe ...`` (drives
    ``latest_release_tag``); ``rev_list`` for ``git rev-list --count ...``
    (the commits-after-tag count). Any other invocation raises, so a test that
    triggers an unexpected git call fails loudly instead of silently passing.
    """

    def _git(*args):
        if args and args[0] == "describe":
            return describe
        if args and args[0] == "rev-list":
            return rev_list
        raise AssertionError(f"unexpected git call: {args!r}")

    return _git


# --- check_dev: the stuck-release invariant -------------------------------


def test_dev_flags_stuck_release(monkeypatch, capsys):
    """version == latest tag AND commits landed after it -> failure (exit 1).

    This is the exact regression the guard exists for: master sitting on an
    already-released version while development continues.
    """
    monkeypatch.setattr(cv, "git", fake_git(describe="v0.9.13", rev_list="276"))
    assert cv.check_dev("0.9.13") == 1
    assert "already released" in capsys.readouterr().err


def test_dev_flags_bare_release_ahead_of_tag(monkeypatch, capsys):
    """A plain release version *ahead* of the latest tag, with commits after it,
    is still a stuck release -> failure (exit 1).

    e.g. the file was bumped to a plain ``0.9.15`` (no ``.dev0``) after tag
    ``v0.9.13`` and development continued. The version differs from the tag, so
    the message names it a plain release rather than an already-released tag.
    """
    monkeypatch.setattr(cv, "git", fake_git(describe="v0.9.13", rev_list="7"))
    assert cv.check_dev("0.9.15") == 1
    assert "a plain release version" in capsys.readouterr().err


def test_dev_allows_transient_release_commit(monkeypatch):
    """version == tag with zero commits after -> allowed (the release commit).

    Keeps the check independent of whether the branch or the tag is pushed
    first; the tagged commit itself is legitimately at the released version.
    """
    monkeypatch.setattr(cv, "git", fake_git(describe="v0.9.14", rev_list="0"))
    assert cv.check_dev("0.9.14") == 0


def test_dev_allows_when_commit_count_unavailable(monkeypatch):
    """git rev-list yielding '' (git unavailable) is treated like zero commits,
    so even a plain release version is allowed rather than falsely flagged.

    Defends the empty-string arm of the ``commits_after in ("", "0")`` guard:
    a plain ``0.9.14`` here must not trip the stuck-release failure.
    """
    monkeypatch.setattr(cv, "git", fake_git(describe="v0.9.14", rev_list=""))
    assert cv.check_dev("0.9.14") == 0


def test_dev_allows_dev_version_with_commits(monkeypatch):
    """A development version with commits after the latest tag passes untouched.

    This is the normal, healthy master state: tag ``v0.9.13`` with the file
    already bumped to ``0.9.14.dev0`` and further commits landed on top.
    """
    monkeypatch.setattr(cv, "git", fake_git(describe="v0.9.13", rev_list="5"))
    assert cv.check_dev("0.9.14.dev0") == 0


def test_dev_allows_no_release_tag(monkeypatch):
    """With no vX.Y.Z tag reachable, there is nothing to be stuck on -> pass."""
    monkeypatch.setattr(cv, "git", fake_git(describe=""))
    assert cv.check_dev("0.9.14.dev0") == 0


# --- check_release: the tag-publish invariant -----------------------------


@pytest.mark.parametrize(
    "version, tag",
    [
        ("0.9.14", "v0.9.14"),  # plain release matches its tag
        ("0.9.7.b0", "v0.9.7.b0"),  # prerelease is fine when the tag matches
        ("1.0.0rc1", "v1.0.0rc1"),  # rc likewise: the tag is the source of truth
    ],
)
def test_release_accepts_matching_tag(version, tag):
    """Any declared version equal to version_from_tag(tag) is publishable,
    including prereleases -- the refactored check no longer rejects suffixes."""
    assert cv.check_release(version, tag) == 0


@pytest.mark.parametrize(
    "version, tag",
    [
        ("0.9.14.dev0", "v0.9.14"),  # dev suffix no longer matches a plain tag
        ("0.9.14", "v0.9.15"),  # plain version, wrong tag
        ("0.9.13", "v0.9.14"),  # off-by-one point release
    ],
)
def test_release_rejects_tag_mismatch(version, tag, capsys):
    """A declared version that differs from version_from_tag(tag) -> failure."""
    assert cv.check_release(version, tag) == 1
    assert "does not match tag" in capsys.readouterr().err


# --- is_release_version: the plain-release classifier ---------------------


@pytest.mark.parametrize(
    "version, expected",
    [
        ("0.9.14", True),  # plain final release
        ("1.0.0", True),  # plain final release
        ("0.9.14.dev0", False),  # post-release development
        ("1.0.0rc1", False),  # release candidate
        ("0.9.7.b0", False),  # beta pre-release
        ("1.0.0+abc", False),  # local version identifier
        ("1.0c1", False),  # only two numeric components
    ],
)
def test_is_release_version(version, expected):
    assert cv.is_release_version(version) is expected


# --- version_from_tag: strip the leading tag 'v' --------------------------


@pytest.mark.parametrize(
    "tag, expected",
    [
        ("v0.9.14", "0.9.14"),  # leading v stripped
        ("v0.9.7.b0", "0.9.7.b0"),  # only the leading v; suffix preserved
        ("0.9.14", "0.9.14"),  # no prefix -> unchanged
    ],
)
def test_version_from_tag(tag, expected):
    assert cv.version_from_tag(tag) == expected


# --- main(): mode dispatch + declared-version wiring ----------------------


def test_main_dev_mode_reports_stuck_release(monkeypatch):
    """Default (dev) mode reads the declared version and detects a stuck one."""
    monkeypatch.setattr(cv, "read_declared_version", lambda: "0.9.13")
    monkeypatch.setattr(cv, "git", fake_git(describe="v0.9.13", rev_list="5"))
    assert cv.main([]) == 1


def test_main_release_mode_matches_tag(monkeypatch):
    """Release mode compares the declared version against --tag and passes."""
    monkeypatch.setattr(cv, "read_declared_version", lambda: "0.9.14")
    assert cv.main(["--mode", "release", "--tag", "v0.9.14"]) == 0


def test_main_release_mode_reports_mismatch(monkeypatch):
    """Release mode wires --tag into check_release: a mismatch fails (exit 1)."""
    monkeypatch.setattr(cv, "read_declared_version", lambda: "0.9.14")
    assert cv.main(["--mode", "release", "--tag", "v0.9.15"]) == 1


def test_main_release_mode_requires_tag(monkeypatch):
    """Release mode without --tag is a usage error (argparse exits nonzero)."""
    monkeypatch.setattr(cv, "read_declared_version", lambda: "0.9.14")
    with pytest.raises(SystemExit) as exc:
        cv.main(["--mode", "release"])
    assert exc.value.code != 0


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-v"]))
