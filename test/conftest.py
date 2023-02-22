import pytest

# Change pytest working directory to test case directory
# https://stackoverflow.com/a/62055409
@pytest.fixture(autouse=True)
def change_test_dir(request, monkeypatch):
    monkeypatch.chdir(request.fspath.dirname)
