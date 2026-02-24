# -*- coding: utf-8 -*-
"""
Unit tests for radis.api.exomolapi functions.
Tests pure parsing functions using synthetic temp files — no network, no GPU.
"""
import numpy as np
import pandas as pd
import pytest


@pytest.mark.fast
class TestMolnameConversion:
    """Tests for exact_molname_exomol_to_simple_molname and helpers."""

    def test_simple_molecule(self):
        from radis.api.exomolapi import exact_molname_exomol_to_simple_molname

        assert exact_molname_exomol_to_simple_molname("12C-1H4") == "CH4"

    def test_multi_atom(self):
        from radis.api.exomolapi import exact_molname_exomol_to_simple_molname

        assert exact_molname_exomol_to_simple_molname("23Na-16O-1H") == "NaOH"

    def test_special_name(self):
        from radis.api.exomolapi import exact_molname_exomol_to_simple_molname

        assert exact_molname_exomol_to_simple_molname("HeH_p") == "HeH_p"

    def test_hho_to_h2o(self):
        from radis.api.exomolapi import exact_molname_exomol_to_simple_molname

        assert exact_molname_exomol_to_simple_molname("1H-2H-16O") == "H2O"

    def test_unconvertable_falls_back(self, capsys):
        from radis.api.exomolapi import exact_molname_exomol_to_simple_molname

        result = exact_molname_exomol_to_simple_molname("trans-31P2-1H-2H")
        assert result == "trans-31P2-1H-2H"


@pytest.mark.fast
class TestComputeWavenumberRanges:
    """Tests for compute_wavenumber_ranges and wavenumber_tag."""

    def test_single_trans_file(self):
        from radis.api.exomolapi import compute_wavenumber_ranges

        result = compute_wavenumber_ranges(1, 10000.0)
        assert result is None

    def test_multiple_trans_files(self):
        from radis.api.exomolapi import compute_wavenumber_ranges

        result = compute_wavenumber_ranges(4, 20000.0)
        np.testing.assert_array_equal(result, [0, 5000, 10000, 15000, 20000])

    def test_wavenumber_tag_none(self):
        from radis.api.exomolapi import wavenumber_tag

        assert wavenumber_tag(None) == ""

    def test_wavenumber_tag_values(self):
        from radis.api.exomolapi import wavenumber_tag

        numinf = np.array([0, 5000, 10000])
        tags = wavenumber_tag(numinf)
        assert tags == ["00000-05000", "05000-10000"]


@pytest.mark.fast
class TestReadDef:
    """Tests for read_def — parsing ExoMol definition files."""

    def test_basic_def_file(self, tmp_path):
        from radis.api.exomolapi import read_def

        # Create a minimal .def file
        content = """0.07  # Default value of Lorentzian half-width
0.5  # Default value of temperature exponent
1  # No. of transition files
20000.0  # Maximum wavenumber (in cm-1)
12.0 0.00199  # Isotopologue mass (Da) and (kg)
0  # Lifetime availability
0  # Lande g-factor availability
0  # Uncertainty availability
v  # Quantum label
"""
        def_file = tmp_path / "test.def"
        def_file.write_text(content)

        result = read_def(def_file)

        assert result["alpha_ref"] == 0.07
        assert result["n_Texp"] == 0.5
        assert result["molmass"] == 12.0
        assert result["numinf"] is None  # single trans file
        assert result["numtag"] == ""
        assert result["quantum_labels"] == ["v"]
        assert result["lifetime"] is False
        assert result["Landé"] is False
        assert result["unc"] is False

    def test_multi_trans_def_file(self, tmp_path):
        from radis.api.exomolapi import read_def

        content = """0.05  # Default value of Lorentzian half-width
0.4  # Default value of temperature exponent
3  # No. of transition files
15000.0  # Maximum wavenumber (in cm-1)
28.0 0.0046  # Isotopologue mass (Da) and (kg)
1  # Lifetime availability
1  # Lande g-factor availability
1  # Uncertainty availability
v  # Quantum label
"""
        def_file = tmp_path / "test.def"
        def_file.write_text(content)

        result = read_def(def_file)

        assert result["alpha_ref"] == 0.05
        assert result["n_Texp"] == 0.4
        assert result["numinf"] is not None
        assert len(result["numinf"]) == 4  # ntransf + 1
        assert result["lifetime"] is True
        assert result["Landé"] is True
        assert result["unc"] is True


@pytest.mark.fast
class TestReadPf:
    """Tests for read_pf — parsing partition function files."""

    def test_basic_pf_file(self, tmp_path):
        from radis.api.exomolapi import read_pf

        content = "100  1.23\n200  4.56\n300  7.89\n"
        pf_file = tmp_path / "test.pf"
        pf_file.write_text(content)

        dat = read_pf(pf_file)
        assert list(dat.columns) == ["T", "QT"]
        assert len(dat) == 3
        assert dat["T"].iloc[0] == 100
        assert dat["QT"].iloc[2] == pytest.approx(7.89)


@pytest.mark.fast
class TestReadTrans:
    """Tests for read_trans — parsing transition files."""

    def test_csv_engine(self, tmp_path):
        from radis.api.exomolapi import read_trans

        content = "1 2 1.5e-10 1000.5\n3 4 2.0e-8 2000.1\n"
        trans_file = tmp_path / "test.trans"
        trans_file.write_text(content)

        dat = read_trans(trans_file, engine="csv")
        assert len(dat) == 2
        assert list(dat.columns) == ["i_upper", "i_lower", "A", "nu_lines"]
        assert dat["i_upper"].iloc[0] == 1
        assert dat["nu_lines"].iloc[1] == pytest.approx(2000.1)

    def test_invalid_engine_raises(self, tmp_path):
        from radis.api.exomolapi import read_trans

        trans_file = tmp_path / "test.trans"
        trans_file.write_text("1 2 1e-10 1000\n")

        with pytest.raises(NotImplementedError):
            read_trans(trans_file, engine="invalid")


@pytest.mark.fast
class TestReadStates:
    """Tests for read_states — parsing state files."""

    def test_csv_skip_optional(self, tmp_path):
        from radis.api.exomolapi import read_states

        content = "1 0.0 1 0\n2 1.448 3 1\n3 4.345 5 2\n"
        states_file = tmp_path / "test.states"
        states_file.write_text(content)

        dic_def = {
            "quantum_labels": [],
            "unc": False,
            "lifetime": False,
            "Landé": False,
        }

        dat = read_states(states_file, dic_def, engine="csv", skip_optional_data=True)
        assert len(dat) == 3
        assert list(dat.columns) == ["i", "E", "g", "J"]
        assert dat["E"].iloc[1] == pytest.approx(1.448)


@pytest.mark.fast
class TestReadBroad:
    """Tests for read_broad and check_code_level."""

    def test_read_broad(self, tmp_path):
        from radis.api.exomolapi import read_broad

        content = (
            "a0 0.07 0.5 0 1 0 0 0 0 0 0 0 0 0 0\na0 0.06 0.4 1 2 0 0 0 0 0 0 0 0 0 0\n"
        )
        broad_file = tmp_path / "test.broad"
        broad_file.write_text(content)

        bdat = read_broad(broad_file, output="pytables")
        assert len(bdat) == 2
        assert bdat["code"].iloc[0] == "a0"
        assert bdat["alpha_ref"].iloc[0] == pytest.approx(0.07)

    def test_check_code_level_a0(self, tmp_path):
        from radis.api.exomolapi import check_code_level, read_broad

        content = (
            "a0 0.07 0.5 0 1 0 0 0 0 0 0 0 0 0 0\na0 0.06 0.4 1 2 0 0 0 0 0 0 0 0 0 0\n"
        )
        broad_file = tmp_path / "test.broad"
        broad_file.write_text(content)

        bdat = read_broad(broad_file)
        assert check_code_level(bdat) == "a0"

    def test_check_code_level_a1(self):
        from radis.api.exomolapi import check_code_level

        bdat = pd.DataFrame({"code": ["a1", "a1"]})
        assert check_code_level(bdat) == "a1"

    def test_check_code_level_mixed(self):
        from radis.api.exomolapi import check_code_level

        bdat = pd.DataFrame({"code": ["a0", "a1"]})
        assert check_code_level(bdat) == "a1"

    def test_check_code_level_m0(self):
        from radis.api.exomolapi import check_code_level

        bdat = pd.DataFrame({"code": ["m0", "m0"]})
        assert check_code_level(bdat) == "m0"

    def test_check_code_level_unknown(self):
        from radis.api.exomolapi import check_code_level

        bdat = pd.DataFrame({"code": ["x0", "x1"]})
        assert check_code_level(bdat) is None


@pytest.mark.fast
class TestMakeJ2b:
    """Tests for make_j2b broadening parameter mapping."""

    def test_basic_j2b(self):
        from radis.api.exomolapi import make_j2b

        bdat = pd.DataFrame(
            {
                "code": ["a0", "a0", "a0"],
                "alpha_ref": [0.08, 0.06, 0.04],
                "n_Texp": [0.6, 0.5, 0.4],
                "jlower": [0, 1, 2],
                "jupper": [1, 2, 3],
            }
        )
        j2a, j2n = make_j2b(bdat, alpha_ref_default=0.07, n_Texp_default=0.5)

        assert len(j2a) == 3
        assert j2a[0] == pytest.approx(0.08)
        assert j2a[1] == pytest.approx(0.06)
        assert j2n[2] == pytest.approx(0.4)


@pytest.mark.fast
class TestPickupGE:
    """Tests for pickup_gE — mapping states into transitions."""

    def test_pytables_engine(self):
        from radis.api.exomolapi import pickup_gE

        states = pd.DataFrame(
            {
                "i": [1, 2, 3, 4],
                "E": [0.0, 100.0, 200.0, 300.0],
                "g": [1, 3, 5, 7],
                "J": [0, 1, 2, 3],
            }
        )
        trans = pd.DataFrame(
            {
                "i_upper": [2, 3],
                "i_lower": [1, 2],
                "A": [1e-10, 2e-8],
                "nu_lines": [100.0, 100.0],
            }
        )
        dic_def = {"quantum_labels": []}

        result = pickup_gE(
            states, trans, dic_def, skip_optional_data=True, engine="pytables"
        )

        assert "gup" in result.columns
        assert "elower" in result.columns
        assert "jlower" in result.columns
        assert "jupper" in result.columns
        assert result["gup"].iloc[0] == 3  # state i=2 has g=3
        assert result["elower"].iloc[0] == pytest.approx(0.0)  # state i=1 has E=0.0
        assert result["jlower"].iloc[1] == 1  # state i=2 has J=1


@pytest.mark.fast
class TestWavenumberRangeHDO:
    """Test wavenumber_range_HDO_VTT."""

    def test_hdo_vtt_returns_array(self):
        from radis.api.exomolapi import wavenumber_range_HDO_VTT

        result = wavenumber_range_HDO_VTT()
        assert isinstance(result, np.ndarray)
        assert result[0] == 0.0
        assert result[-1] == 26000.0
        assert len(result) == 17


if __name__ == "__main__":
    pytest.main([__file__])
