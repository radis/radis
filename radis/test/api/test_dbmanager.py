import uuid
import bz2
from pathlib import Path

import pytest
import requests

from radis.api.dbmanager import DatabaseManager


class _DummyDatabaseManager(DatabaseManager):
    def __init__(self, local_databases):
        super().__init__(
            name=f"DBMANAGER-TEST-{uuid.uuid4()}",
            molecule="X",
            local_databases=str(local_databases),
            engine="pytables",
            verbose=False,
            parallel=False,
        )
        self.downloadable = False
        self.parsed_payload = None

    def parse_to_local_file(
        self,
        opener,
        urlname,
        local_file,
        pbar_active=True,
        pbar_t0=0,
        pbar_Ntot_estimate_factor=None,
        pbar_Nlines_already=0,
        pbar_last=True,
    ):
        with opener.open(urlname) as stream:
            self.parsed_payload = stream.read()
        Path(local_file).write_bytes(self.parsed_payload)
        return len(self.parsed_payload)


class _FakeResponse:
    def __init__(self, status_code, headers, chunks):
        self.status_code = status_code
        self.headers = headers
        self._chunks = chunks

    def raise_for_status(self):
        if self.status_code >= 400:
            raise requests.HTTPError(f"HTTP {self.status_code}")

    def iter_content(self, chunk_size=8192):
        del chunk_size
        for chunk in self._chunks:
            yield chunk


class _FakeSession:
    def __init__(self, head_response, get_response):
        self._head_response = head_response
        self._get_response = get_response

    def head(self, urlname, headers=None, allow_redirects=True):
        del urlname, headers, allow_redirects
        return self._head_response

    def get(self, urlname, headers=None, stream=True, allow_redirects=True):
        del urlname, headers, stream, allow_redirects
        return self._get_response


def _run_single_download(manager, local_file):
    manager.download_and_parse(
        ["https://example.org/testfile.bin"],
        [str(local_file)],
        N_files_total=1,
    )


@pytest.mark.fast
def test_download_and_parse_allows_html_head_when_get_is_binary(monkeypatch, tmp_path):
    monkeypatch.chdir(tmp_path)
    local_file = tmp_path / "output.bin"
    manager = _DummyDatabaseManager(tmp_path / "db")

    session = _FakeSession(
        head_response=_FakeResponse(
            status_code=200,
            headers={"content-type": "text/html"},
            chunks=[],
        ),
        get_response=_FakeResponse(
            status_code=200,
            headers={"content-type": "application/octet-stream"},
            chunks=[b"radis-test-binary-content"],
        ),
    )
    monkeypatch.setattr("radis.api.dbmanager.requests.Session", lambda: session)

    _run_single_download(manager, local_file)

    assert manager.parsed_payload == b"radis-test-binary-content"
    assert local_file.read_bytes() == b"radis-test-binary-content"


@pytest.mark.fast
def test_download_and_parse_raises_for_html_get_page(monkeypatch, tmp_path):
    monkeypatch.chdir(tmp_path)
    local_file = tmp_path / "output.bin"
    manager = _DummyDatabaseManager(tmp_path / "db")

    session = _FakeSession(
        head_response=_FakeResponse(
            status_code=200,
            headers={"content-type": "application/octet-stream"},
            chunks=[],
        ),
        get_response=_FakeResponse(
            status_code=200,
            headers={"content-type": "text/html"},
            chunks=[b"<!DOCTYPE html><html><body>login</body></html>"],
        ),
    )
    monkeypatch.setattr("radis.api.dbmanager.requests.Session", lambda: session)

    with pytest.raises(requests.HTTPError, match="Received an HTML page"):
        _run_single_download(manager, local_file)


@pytest.mark.fast
def test_download_and_parse_raises_for_head_http_error(monkeypatch, tmp_path):
    monkeypatch.chdir(tmp_path)
    local_file = tmp_path / "output.bin"
    manager = _DummyDatabaseManager(tmp_path / "db")

    class _NoGetSession:
        def head(self, urlname, headers=None, allow_redirects=True):
            del urlname, headers, allow_redirects
            return _FakeResponse(status_code=403, headers={}, chunks=[])

        def get(self, urlname, headers=None, stream=True, allow_redirects=True):
            del urlname, headers, stream, allow_redirects
            raise AssertionError("GET must not be called when HEAD fails")

    monkeypatch.setattr("radis.api.dbmanager.requests.Session", _NoGetSession)

    with pytest.raises(requests.HTTPError, match="status code 403"):
        _run_single_download(manager, local_file)


@pytest.mark.fast
def test_download_and_parse_allows_redirect_like_head_status(monkeypatch, tmp_path):
    monkeypatch.chdir(tmp_path)
    local_file = tmp_path / "output.bin"
    manager = _DummyDatabaseManager(tmp_path / "db")

    session = _FakeSession(
        head_response=_FakeResponse(status_code=302, headers={"content-type": "text/html"}, chunks=[]),
        get_response=_FakeResponse(
            status_code=200,
            headers={"content-type": "application/octet-stream"},
            chunks=[b"redirect-followed-binary"],
        ),
    )
    monkeypatch.setattr("radis.api.dbmanager.requests.Session", lambda: session)

    _run_single_download(manager, local_file)
    assert manager.parsed_payload == b"redirect-followed-binary"


@pytest.mark.fast
def test_download_and_parse_uses_hitran_session_for_hitemp_urls(monkeypatch, tmp_path):
    monkeypatch.chdir(tmp_path)
    local_file = tmp_path / "output.bin"
    manager = _DummyDatabaseManager(tmp_path / "db")

    session = _FakeSession(
        head_response=_FakeResponse(
            status_code=200,
            headers={"content-type": "text/html"},
            chunks=[],
        ),
        get_response=_FakeResponse(
            status_code=200,
            headers={"content-type": "application/octet-stream"},
            chunks=[bz2.compress(b"hitemp-download-content")],
        ),
    )
    calls = {"login": 0}

    def _fake_login_to_hitran(verbose=False):
        del verbose
        calls["login"] += 1
        return session

    def _unexpected_requests_session():
        raise AssertionError("requests.Session() should not be used for HITEMP URLs")

    monkeypatch.setattr("radis.api.hitempapi.login_to_hitran", _fake_login_to_hitran)
    monkeypatch.setattr("radis.api.dbmanager.requests.Session", _unexpected_requests_session)

    manager.download_and_parse(
        ["https://hitran.org/files/HITEMP/bzip2format/06_HITEMP2020.par.bz2"],
        [str(local_file)],
        N_files_total=1,
    )

    assert calls["login"] == 1
    assert manager.parsed_payload == b"hitemp-download-content"
