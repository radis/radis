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


@pytest.fixture
def manager(tmp_path):
    return _DummyDownloadManager(tmp_path)


@pytest.fixture
def run_download(tmp_path, monkeypatch):
    def _run(
        manager, session, *, urlname="https://example.org/test.par.bz2", hitemp=False
    ):
        if hitemp:
            import radis.api.hitempapi as hitempapi

            monkeypatch.setattr(
                hitempapi, "login_to_hitran", lambda verbose=False: session
            )
        else:
            monkeypatch.setattr(requests, "Session", lambda: session)

        monkeypatch.chdir(tmp_path)
        manager.download_and_parse(
            urlnames=[urlname],
            local_files=[str(tmp_path / "dummy_output.hdf5")],
            N_files_total=1,
        )

    return _run


@pytest.mark.fast
def test_hitemp_head_html_get_binary_does_not_fail(manager, run_download):
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
    run_download(
        manager,
        session,
        urlname="https://hitran.org/files/HITEMP/bzip2format/06_HITEMP2020.par.bz2",
        hitemp=True,
    )

    assert manager.parse_calls == 1
    assert len(session.head_calls) == 1
    assert len(session.get_calls) == 1


@pytest.mark.fast
@pytest.mark.parametrize(
    "content_type,chunks",
    [
        ("text/html", [b"<html><body>login</body></html>"]),
        ("text/html; charset=utf-8", [b"\x00\x00"]),
        ("application/octet-stream", [b"<!doctype html><html>login</html>"]),
    ],
)
def test_get_invalid_payload_raises_httperror(
    manager, run_download, content_type, chunks
):
    session = _FakeSession(
        head_response=_FakeResponse(
            status_code=200, headers={"content-type": "text/html"}
        ),
        get_response=_FakeResponse(
            status_code=200,
            headers={"content-type": content_type, "content-length": "32"},
            chunks=chunks,
        ),
    )

    with pytest.raises(requests.HTTPError, match="HTML content from GET request"):
        run_download(manager, session)

    assert manager.parse_calls == 0


@pytest.mark.fast
@pytest.mark.parametrize("status_code", [404, 302])
def test_head_non_200_raises_before_get(manager, run_download, status_code):
    session = _FakeSession(
        head_response=_FakeResponse(
            status_code=status_code, headers={"content-type": "text/html"}
        ),
        get_response=_FakeResponse(
            status_code=200,
            headers={"content-type": "application/octet-stream"},
            chunks=[b"data"],
        ),
    )

    with pytest.raises(requests.HTTPError, match="HEAD request"):
        run_download(manager, session)

    assert manager.parse_calls == 0
    assert len(session.get_calls) == 0
