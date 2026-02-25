# -*- coding: utf-8 -*-
import pytest
import requests

from radis.api.dbmanager import DatabaseManager


class _FakeResponse:
    def __init__(self, status_code=200, headers=None, chunks=None):
        self.status_code = status_code
        self.headers = headers or {}
        self._chunks = chunks or []

    def raise_for_status(self):
        if self.status_code >= 400:
            raise requests.HTTPError(f"status {self.status_code}")

    def iter_content(self, chunk_size=8192):
        for chunk in self._chunks:
            yield chunk


class _FakeSession:
    def __init__(self, head_response, get_response):
        self._head_response = head_response
        self._get_response = get_response
        self.head_calls = []
        self.get_calls = []

    def head(self, *args, **kwargs):
        self.head_calls.append((args, kwargs))
        return self._head_response

    def get(self, *args, **kwargs):
        self.get_calls.append((args, kwargs))
        return self._get_response


class _DummyDownloadManager(DatabaseManager):
    def __init__(self, local_databases):
        super().__init__(
            name="TEST-DOWNLOAD",
            molecule="CH4",
            local_databases=str(local_databases),
            engine="pytables",
            verbose=False,
            parallel=False,
        )
        self.parse_calls = 0

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
        self.parse_calls += 1
        return 1


@pytest.mark.fast
def test_hitemp_head_html_get_binary_does_not_fail(tmp_path, monkeypatch):
    manager = _DummyDownloadManager(tmp_path)
    session = _FakeSession(
        head_response=_FakeResponse(
            status_code=200,
            headers={"content-type": "text/html"},
        ),
        get_response=_FakeResponse(
            status_code=200,
            headers={"content-type": "application/octet-stream", "content-length": "8"},
            chunks=[b"BZh9", b"DATA"],
        ),
    )

    import radis.api.hitempapi as hitempapi

    monkeypatch.setattr(hitempapi, "login_to_hitran", lambda verbose=False: session)
    monkeypatch.chdir(tmp_path)

    manager.download_and_parse(
        urlnames=["https://hitran.org/files/HITEMP/bzip2format/06_HITEMP2020.par.bz2"],
        local_files=[str(tmp_path / "dummy_output.hdf5")],
        N_files_total=1,
    )

    assert manager.parse_calls == 1
    assert len(session.head_calls) == 1
    assert len(session.get_calls) == 1


@pytest.mark.fast
def test_get_html_payload_raises_httperror(tmp_path, monkeypatch):
    manager = _DummyDownloadManager(tmp_path)
    session = _FakeSession(
        head_response=_FakeResponse(
            status_code=200, headers={"content-type": "text/html"}
        ),
        get_response=_FakeResponse(
            status_code=200,
            headers={"content-type": "text/html", "content-length": "32"},
            chunks=[b"<html><body>login</body></html>"],
        ),
    )

    monkeypatch.setattr(requests, "Session", lambda: session)
    monkeypatch.chdir(tmp_path)

    with pytest.raises(requests.HTTPError, match="HTML content from GET request"):
        manager.download_and_parse(
            urlnames=["https://example.org/test.par.bz2"],
            local_files=[str(tmp_path / "dummy_output.hdf5")],
            N_files_total=1,
        )

    assert manager.parse_calls == 0


@pytest.mark.fast
def test_head_http_error_raises_before_get(tmp_path, monkeypatch):
    manager = _DummyDownloadManager(tmp_path)
    session = _FakeSession(
        head_response=_FakeResponse(
            status_code=404, headers={"content-type": "text/html"}
        ),
        get_response=_FakeResponse(
            status_code=200,
            headers={"content-type": "application/octet-stream"},
            chunks=[b"data"],
        ),
    )

    monkeypatch.setattr(requests, "Session", lambda: session)

    with pytest.raises(requests.HTTPError, match="HEAD request"):
        manager.download_and_parse(
            urlnames=["https://example.org/test.par.bz2"],
            local_files=[str(tmp_path / "dummy_output.hdf5")],
            N_files_total=1,
        )

    assert manager.parse_calls == 0
    assert len(session.get_calls) == 0


@pytest.mark.fast
def test_head_redirect_status_is_accepted(tmp_path, monkeypatch):
    manager = _DummyDownloadManager(tmp_path)
    session = _FakeSession(
        head_response=_FakeResponse(
            status_code=302,
            headers={
                "content-type": "text/html",
                "location": "https://example.org/final",
            },
        ),
        get_response=_FakeResponse(
            status_code=200,
            headers={"content-type": "application/octet-stream", "content-length": "4"},
            chunks=[b"data"],
        ),
    )

    monkeypatch.setattr(requests, "Session", lambda: session)
    monkeypatch.chdir(tmp_path)

    manager.download_and_parse(
        urlnames=["https://example.org/test.par.bz2"],
        local_files=[str(tmp_path / "dummy_output.hdf5")],
        N_files_total=1,
    )

    assert manager.parse_calls == 1
