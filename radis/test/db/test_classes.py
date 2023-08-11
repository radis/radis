from radis.db.classes import is_atom
def test_is_atom():
    
    assert is_atom("Fe_I")
    assert is_atom("H_II")
    assert not is_atom("CO2")

