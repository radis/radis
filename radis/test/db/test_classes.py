from radis.db.classes import is_atom,to_conventional_name
def test_is_atom():
    
    assert is_atom("Fe_I")
    assert is_atom("H_II")
    assert not is_atom("CO2")

def test_to_conventional_name():
    assert to_conventional_name("Fe")=="Fe_I"
    assert to_conventional_name("H+")=="H_II"
    assert to_conventional_name("CO2")=="CO2"

