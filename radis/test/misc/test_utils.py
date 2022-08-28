from os.path import basename, dirname, join

from radis.misc.utils import get_files_from_regex
from radis.test.utils import getTestFile

# Test to compare basenames of files whose names match the regular expression passed.
# TEST_FOLDER_PATH = dirname(abspath(dirname(__file__)))
TEST_FOLDER_PATH = dirname(getTestFile("."))


def test_files_from_regex_1(*args, **kwargs):
    """
    Check if test_file_from_regex() is returning a list of correct path values
    This test case uses '*.par' as the regular expression
    path1 generates a path-like string for '*.par'
    path_test contains all the files matching '*.par' regular expression
    path_actual contains the expected list of file names which should be returned by get_files_from_regex
    Lists are looped through to check if every element in path_test and path_actual match, if not an AssertError is thrown
    """

    path = join(TEST_FOLDER_PATH, "*.par")
    # path_test is for finding all the files with .par extension having regular expression *.par
    path_test = get_files_from_regex(path)

    # This list will have to be updated if more .par files are added/removed in/from the folder
    # This list is what is the expected output from get_files_from_regex()
    path_actual = [
        "geisa_CO_fragment.par",
        "geisa_CO2_fragment.par",
        "geisa_H2O_fragment.par",
        "geisa_O2_fragment.par",
        "hitran_co2_626_bandhead_4165_4200nm.par",
        "hitran_2016_H2O_2iso_2000_2100cm.par",
        "hitran_CO_fragment.par",
        "hitran_co_3iso_2000_2300cm.par",
        "hitran_CO2_fragment.par",
    ]

    assert len(path_test) == len(path_actual)

    path_actual.sort()
    path_test.sort()

    for i in range(
        len(path_actual)
    ):  # finding basename because the function `get_file_from_regex` returns the absolute path
        path_test[i] = basename(path_test[i])
        assert path_test[i] == path_actual[i]


def test_files_from_regex_2(*args, **kwargs):
    """
    Check if test_file_from_regex() is returning a list of correct path values
    This test case uses '*co2*.par' as the regular expression
    path generates a path-like string for '*co2*.par'
    path_test contains all the files matching '*co2*.par' regular expression
    path_actual contains the expected list of file names which should be returned by get_files_from_regex
    Lists are looped through to check if every element in path_test and path_actual match, if not an AssertError is thrown
    """

    path = join(TEST_FOLDER_PATH, "*co2*.par")

    path_test = get_files_from_regex(path)

    # This list will have to be updated if more files matching regular expression cdsd_* are added/removed in/from the folder
    # This list is what is the expected output from get_files_from_regex()
    path_actual = [
        "geisa_CO2_fragment.par",
        "hitran_CO2_fragment.par",
        "hitran_co2_626_bandhead_4165_4200nm.par",
    ]
    path_actual.sort()
    path_test.sort()

    for i in range(
        len(path_actual)
    ):  # finding basename because the function `get_file_from_regex` returns the absolute path
        path_test[i] = basename(path_test[i])
        assert path_test[i] == path_actual[i]

    # print(path_test2)


def _run_testcases(verbose=True, *args, **kwargs):
    test_files_from_regex_1()
    test_files_from_regex_2()
    return True


if __name__ == "__main__":
    print("test util files:", _run_testcases())
