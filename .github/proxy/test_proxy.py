# import os
import sys
import requests
import urllib3

expected_outcome = [
    [("MISS" if sys.argv[-1] == "pre-build" else "HIT"), "HIT"],
    ["HIT", "HIT"],
    [None, None],
]

headers = {
    "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36"
}

url1 = "https://www.hitran.org/data/Aerosols/Aerosols-2024/Aerosol_Readme_2024.pdf"
url2 = "https://hitran.org/data/Aerosols/Aerosols-2024/Aerosol_Readme_2024.pdf"
url3 = "https://geisa.aeris-data.fr/wp-content-aeris/uploads/sites/33/2021/07/table_abundance_GEISA_2020.pdf"


def format_header(s):
    if s is None:
        return "NONE"
    elif s == "HIT":
        return "HIT "
    else:
        return s


for i, url in enumerate([url1, url2, url3]):

    print(i + 1, "requests", end=" ")
    r = requests.get(url)
    print(r.status_code, end=" ")
    print(format_header(r.headers.get("X-Cache-Status")), end=" ")
    exp_str = format_header(expected_outcome[i][0])
    if (r.headers.get("X-Cache-Status") != expected_outcome[i][0]) or (
        r.status_code != 200
    ):
        print(f"!= {exp_str} (FAIL)")
        sys.exit(1)
    else:
        print(f"== {exp_str} (PASS)")

    print(i + 1, "urllib3 ", end=" ")
    http = urllib3.PoolManager()
    r = http.request("GET", url, headers=headers)
    print(r.status, end=" ")
    print(format_header(r.headers.get("X-Cache-Status")), end=" ")
    exp_str = format_header(expected_outcome[i][1])

    if (r.headers.get("X-Cache-Status") != expected_outcome[i][1]) or (r.status != 200):
        print(f"!= {exp_str} (FAIL)")
        sys.exit(1)
    else:
        print(f"== {exp_str} (PASS)")
