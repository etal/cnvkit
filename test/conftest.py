import pytest


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
