import ast
import inspect

import matplotlib
import pytest

# Use a headless backend for all plotting tests (set once, before any test
# module imports pyplot via cnvlib.scatter/plots).
matplotlib.use("Agg")


# Change pytest working directory to test case directory
# https://stackoverflow.com/a/62055409
@pytest.fixture(autouse=True)
def change_test_dir(request, monkeypatch):
    monkeypatch.chdir(request.fspath.dirname)


def linecount(filename):
    i = -1
    with open(filename) as handle:
        for i, _line in enumerate(handle):
            pass
        return i + 1


def ast_calls_to(func, attr, obj=None):
    """Return the ``ast.Call`` nodes in `func`'s source that call `attr`.

    `obj`, when given, additionally requires the call to be made on that name,
    so ``ast_calls_to(batch.batch_run_sample, "do_fix", "fix")`` matches only
    ``fix.do_fix(...)``. Used by the argument-propagation guards, which assert
    that a high-level entry point forwards an option to the function that
    implements it -- a regression class that otherwise needs a BAM fixture and
    a full pipeline run to detect.

    Only ``ast.Attribute`` calls are matched, so a call reached through an
    aliased import (``from cnvlib.call import do_call as _dc``) is invisible
    here. A guard over a function that calls `attr` more than once must
    therefore pin the number of invocations it expects: an existence check
    would still pass on whichever sites remain visible, leaving the aliased
    one unexamined.
    """
    calls = []
    for node in ast.walk(ast.parse(inspect.getsource(func))):
        if not isinstance(node, ast.Call) or not isinstance(node.func, ast.Attribute):
            continue
        if node.func.attr != attr:
            continue
        if obj is not None and not (
            isinstance(node.func.value, ast.Name) and node.func.value.id == obj
        ):
            continue
        calls.append(node)
    return calls


def ast_submit_calls(func, target):
    """Return the ``pool.submit(target, ...)`` call nodes in `func`'s source.

    `target` is matched by the name written at the call site, so the caveat on
    `ast_calls_to` applies here too: submitting the function under an alias
    hides the call, and a guard over more than one submit site must pin the
    count it expects.
    """
    calls = []
    for call in ast_calls_to(func, "submit"):
        if not call.args:
            continue
        first = call.args[0]
        name = first.id if isinstance(first, ast.Name) else getattr(first, "attr", None)
        if name == target:
            calls.append(call)
    return calls


def call_arg_names(call):
    """Names, attribute names and keywords supplied at a call site.

    A propagation guard can then accept ``f(no_overlap)``, ``f(args.no_overlap)``
    or ``f(no_overlap=...)`` interchangeably.
    """
    return (
        {a.id for a in call.args if isinstance(a, ast.Name)}
        | {a.attr for a in call.args if isinstance(a, ast.Attribute)}
        | {kw.arg for kw in call.keywords}
    )
