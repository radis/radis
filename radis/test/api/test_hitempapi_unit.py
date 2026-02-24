# -*- coding: utf-8 -*-
"""
Unit tests for radis.api.hitempapi functions.
Tests pure functions using synthetic data — no network, no GPU.
"""
import json
import os

import numpy as np
import pytest


@pytest.mark.fast
class TestKeepOnlyRelevant:
    """Tests for keep_only_relevant — filtering files by wavenumber range."""

    def test_both_bounds(self):
        from radis.api.hitempapi import keep_only_relevant

        files = [
            "/data/CO2-02_00000-00250_HITEMP.h5",
            "/data/CO2-02_00250-00500_HITEMP.h5",
            "/data/CO2-02_00500-00750_HITEMP.h5",
            "/data/CO2-02_00750-01000_HITEMP.h5",
        ]
        relevant, wmin, wmax = keep_only_relevant(
            files, wavenum_min=400, wavenum_max=600, verbose=False
        )
        assert len(relevant) == 2
        assert "/data/CO2-02_00250-00500_HITEMP.h5" in relevant
        assert "/data/CO2-02_00500-00750_HITEMP.h5" in relevant
        assert wmin == 250
        assert wmax == 750

    def test_no_bounds(self):
        from radis.api.hitempapi import keep_only_relevant

        files = [
            "/data/CO2-02_00000-00250_HITEMP.h5",
            "/data/CO2-02_00250-00500_HITEMP.h5",
        ]
        relevant, _, _ = keep_only_relevant(files, verbose=False)
        assert len(relevant) == 2

    def test_only_min_bound(self):
        from radis.api.hitempapi import keep_only_relevant

        files = [
            "/data/CO2-02_00000-00250_HITEMP.h5",
            "/data/CO2-02_00250-00500_HITEMP.h5",
            "/data/CO2-02_00500-00750_HITEMP.h5",
        ]
        relevant, _, _ = keep_only_relevant(files, wavenum_min=400, verbose=False)
        assert len(relevant) == 2  # 250-500 (wmax=500 > 400) and 500-750

    def test_only_max_bound(self):
        from radis.api.hitempapi import keep_only_relevant

        files = [
            "/data/CO2-02_00000-00250_HITEMP.h5",
            "/data/CO2-02_00250-00500_HITEMP.h5",
            "/data/CO2-02_00500-00750_HITEMP.h5",
        ]
        relevant, _, _ = keep_only_relevant(files, wavenum_max=300, verbose=False)
        assert len(relevant) == 2  # 0-250 (wmin=0 < 300) and 250-500 (wmin=250 < 300)

    def test_no_relevant_files(self):
        from radis.api.hitempapi import keep_only_relevant

        files = [
            "/data/CO2-02_00000-00250_HITEMP.h5",
        ]
        relevant, wmin, wmax = keep_only_relevant(
            files, wavenum_min=500, wavenum_max=600, verbose=False
        )
        assert len(relevant) == 0

    def test_verbose_single_file(self, capsys):
        from radis.api.hitempapi import keep_only_relevant

        files = ["/data/CO2-02_00000-00250_HITEMP.h5"]
        relevant, _, _ = keep_only_relevant(
            files, wavenum_min=100, wavenum_max=200, verbose=True
        )
        captured = capsys.readouterr()
        assert "Keep only relevant input file" in captured.out


@pytest.mark.fast
class TestGetLast:
    """Tests for get_last — byte array parsing."""

    def test_basic(self):
        from radis.api.hitempapi import get_last

        # Array with non-empty items followed by empty ones
        b = np.array(["hello world", "another line", "", ""], dtype=object)
        result = get_last(b)
        assert len(result) == 2
        assert result[0] == "hello world"


@pytest.mark.fast
class TestRunningInSpyder:
    """Tests for running_in_spyder."""

    def test_not_in_spyder(self):
        from radis.api.hitempapi import running_in_spyder

        # Ensure SPYDER_ARGS is not set
        os.environ.pop("SPYDER_ARGS", None)
        assert running_in_spyder() is False


@pytest.mark.fast
class TestReadConfig:
    """Tests for read_config."""

    def test_read_existing_config(self, tmp_path, monkeypatch):
        from radis.api.hitempapi import read_config

        config_file = tmp_path / "radis.json"
        config_file.write_text(json.dumps({"key": "value"}))
        monkeypatch.setattr("radis.api.hitempapi.CONFIG_PATH_JSON", str(config_file))

        result = read_config()
        assert result == {"key": "value"}

    def test_read_nonexistent_config(self, monkeypatch):
        from radis.api.hitempapi import read_config

        monkeypatch.setattr(
            "radis.api.hitempapi.CONFIG_PATH_JSON", "/nonexistent/path.json"
        )
        result = read_config()
        assert result == {}


@pytest.mark.fast
class TestEncryptionRoundtrip:
    """Tests for encrypt_password / decrypt_password roundtrip."""

    def test_roundtrip(self, tmp_path, monkeypatch):
        from radis.api.hitempapi import decrypt_password, encrypt_password

        # Use a temporary config file so we don't modify the real one
        config_file = tmp_path / "radis.json"
        config_file.write_text("{}")
        monkeypatch.setattr("radis.api.hitempapi.CONFIG_PATH_JSON", str(config_file))

        original = "my_secret_password_123!"
        encrypted = encrypt_password(original)
        assert encrypted != original
        decrypted = decrypt_password(encrypted)
        assert decrypted == original


@pytest.mark.fast
class TestSetupCredentials:
    """Tests for setup_credentials."""

    def test_env_vars_set(self, monkeypatch):
        from radis.api.hitempapi import setup_credentials

        monkeypatch.setenv("HITRAN_EMAIL", "test@example.com")
        monkeypatch.setenv("HITRAN_PASSWORD", "secret123")
        # Remove CI env vars
        monkeypatch.delenv("READTHEDOCS", raising=False)
        monkeypatch.delenv("TRAVIS", raising=False)
        monkeypatch.delenv("GITHUB_ACTIONS", raising=False)

        email, password = setup_credentials()
        assert email == "test@example.com"
        assert password == "secret123"

    def test_ci_env_no_email_warns(self, monkeypatch):
        from radis.api.hitempapi import setup_credentials

        monkeypatch.setenv("GITHUB_ACTIONS", "true")
        monkeypatch.delenv("HITRAN_EMAIL", raising=False)
        monkeypatch.setenv("HITRAN_PASSWORD", "secret123")

        with pytest.warns(UserWarning, match="HITRAN_EMAIL"):
            email, password = setup_credentials()

    def test_ci_env_no_password_warns(self, monkeypatch):
        from radis.api.hitempapi import setup_credentials

        monkeypatch.setenv("TRAVIS", "true")
        monkeypatch.setenv("HITRAN_EMAIL", "test@example.com")
        monkeypatch.delenv("HITRAN_PASSWORD", raising=False)

        with pytest.warns(UserWarning, match="HITRAN_PASSWORD"):
            email, password = setup_credentials()

    def test_pytest_env_no_creds_raises(self, monkeypatch):
        from radis.api.hitempapi import setup_credentials

        monkeypatch.setenv("PYTEST_CURRENT_TEST", "test_something")
        monkeypatch.delenv("HITRAN_EMAIL", raising=False)
        monkeypatch.delenv("HITRAN_PASSWORD", raising=False)
        monkeypatch.delenv("READTHEDOCS", raising=False)
        monkeypatch.delenv("TRAVIS", raising=False)
        monkeypatch.delenv("GITHUB_ACTIONS", raising=False)

        with pytest.raises(OSError, match="HITRAN_EMAIL"):
            setup_credentials()


@pytest.mark.fast
class TestStoreCredentials:
    """Tests for store_credentials."""

    def test_store_and_verify(self, tmp_path, monkeypatch):
        from radis.api.hitempapi import store_credentials

        config_file = tmp_path / "radis.json"
        config_file.write_text("{}")
        monkeypatch.setattr("radis.api.hitempapi.CONFIG_PATH_JSON", str(config_file))

        store_credentials("user@example.com", "password123")

        with open(config_file) as f:
            config = json.load(f)

        assert "credentials" in config
        assert "HITRAN_EMAIL" in config["credentials"]
        assert "HITRAN_PASSWORD" in config["credentials"]
        # Verify they're encrypted (not plaintext)
        assert config["credentials"]["HITRAN_EMAIL"] != "user@example.com"
        assert config["credentials"]["HITRAN_PASSWORD"] != "password123"


@pytest.mark.fast
class TestFcacheFileName:
    """Tests for _fcache_file_name."""

    def test_pytables_engine(self):
        from radis.api.hitempapi import _fcache_file_name

        result = _fcache_file_name("/data/test.par", "pytables")
        assert str(result).endswith(".h5")

    def test_vaex_engine(self):
        from radis.api.hitempapi import _fcache_file_name

        result = _fcache_file_name("/data/test.par", "vaex")
        assert str(result).endswith(".hdf5")


@pytest.mark.fast
class TestLoadCacheFile:
    """Tests for _load_cache_file."""

    def test_nonexistent_returns_none(self):
        from radis.api.hitempapi import _load_cache_file

        result = _load_cache_file("/nonexistent/cache.h5")
        assert result is None


@pytest.mark.fast
class TestGetEncryptionKey:
    """Tests for get_encryption_key."""

    def test_creates_key_if_missing(self, tmp_path, monkeypatch):
        from radis.api.hitempapi import get_encryption_key

        config_file = tmp_path / "radis.json"
        config_file.write_text("{}")
        monkeypatch.setattr("radis.api.hitempapi.CONFIG_PATH_JSON", str(config_file))

        key = get_encryption_key()
        assert isinstance(key, bytes)
        assert len(key) > 0

        # Verify key is stored
        with open(config_file) as f:
            config = json.load(f)
        assert "credentials" in config
        assert "ENCRYPTION_KEY" in config["credentials"]

    def test_returns_existing_key(self, tmp_path, monkeypatch):
        from cryptography.fernet import Fernet

        from radis.api.hitempapi import get_encryption_key

        config_file = tmp_path / "radis.json"
        existing_key = Fernet.generate_key().decode()
        config_file.write_text(
            json.dumps({"credentials": {"ENCRYPTION_KEY": existing_key}})
        )
        monkeypatch.setattr("radis.api.hitempapi.CONFIG_PATH_JSON", str(config_file))

        key = get_encryption_key()
        assert key == existing_key.encode()


if __name__ == "__main__":
    pytest.main([__file__])
