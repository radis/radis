import os

from radis.api.hitempapi import login_to_hitran
from radis.io.hitemp import fetch_hitemp
from radis.test.utils import getTestFile


def test_hitemp_download():
    """Test downloading a small HITEMP database (OH) and verify basic functionality.

    This test:
    1. Tests login to HITRAN
    2. Downloads OH database (smallest HITEMP molecule, ~900kb)
    3. Verifies basic parsing functionality
    """
    print("\n=== HITEMP Download Test ===")
    print("This test requires HITRAN credentials.")

    # Test login
    print("\nTesting HITRAN login...")
    session = login_to_hitran(verbose=True)
    assert session is not None, "Login failed"
    print("Login successful!")

    # Test downloading OH (smallest molecule)
    print("\nDownloading OH database (~900kb)...")
    df = fetch_hitemp(
        "OH",
        local_databases=os.path.join(getTestFile("."), "hitemp"),
        databank_name="HITEMP-OH-TEST",
        verbose=True,
        cache="regen",  # Force redownload
    )
    print(df)
    # Verify basic functionality
    print("\nVerifying downloaded data...")
    assert len(df) == 57019, f"Expected 57019 lines, got {len(df)}"
    assert "wav" in df.columns, "Wavelength column missing"
    assert "int" in df.columns, "Intensity column missing"

    # Test basic wavelength filtering
    df_filtered = fetch_hitemp(
        "OH",
        local_databases=os.path.join(getTestFile("."), "hitemp"),
        databank_name="HITEMP-OH-TEST",
        load_wavenum_min=2000,
        load_wavenum_max=3000,
        verbose=True,
    )
    assert len(df_filtered) < len(df), "Wavelength filtering not working"

    print("\nAll tests passed successfully!")
    print("HITEMP download functionality is working correctly.")


if __name__ == "__main__":
    test_hitemp_download()
