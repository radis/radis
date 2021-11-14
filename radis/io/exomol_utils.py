"""Utils to read ExoMol

Borrowed from the `Exojax <https://github.com/HajimeKawahara/exojax>`__
code (which you should also have a look at !), by @HajimeKawahara, under MIT License.

"""

import re


def e2s(molname_exact):
    """convert the exact molname (used in ExoMol) to the simple molname

    Args:
       molname_exact: the exact molname

    Returns:
       simple molname

    Examples::

       >>> print(e2s("12C-1H4"))
       >>> CH4
       >>> print(e2s("23Na-16O-1H"))
       >>> NaOH
       >>> print(e2s("HeH_p"))
       >>> HeH_p
       >>> print(e2s("trans-31P2-1H-2H")) #not working
       >>> Warning: Exact molname  trans-31P2-1H-2H cannot be converted to simple molname
       >>> trans-31P2-1H-2H

    """

    try:
        t = molname_exact.split("-")
        molname_simple = ""
        for ele in t:
            alp = "".join(re.findall(r"\D", ele))
            num0 = re.split("[A-Z]", ele)[1]
            if num0.isdigit():
                num = num0
            else:
                num = ""
            molname_simple = molname_simple + alp + num
        return molname_simple
    except:
        print(
            "Warning: Exact molname ",
            molname_exact,
            "cannot be converted to simple molname",
        )
        return molname_exact


if __name__ == "__main__":
    print(e2s("12C-1H4"))
    print(e2s("23Na-16O-1H"))
    print(e2s("HeH_p"))
    print(e2s("trans-31P2-1H-2H"))  # not working

    # TODO : move to unitary tests of exomol_utils
    assert e2s("12C-1H4") == "CH4"
    assert e2s("23Na-16O-1H") == "NaOH"
    assert e2s("HeH_p") == "HeH_p"
    assert e2s("trans-31P2-1H-2H") == "trans-31P2-1H-2H"  # convert not working

    assert e2s("12C-16O2") == "CO2"
    assert e2s("13C-16O2") == "CO2"
