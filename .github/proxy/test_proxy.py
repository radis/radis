import os

import requests

print("CA_bundle:", os.environ["REQUESTS_CA_BUNDLE"])
print("${{ github.workspace }} path:", os.path.abspath("."))
print("${{ github.workspace }} contents:", os.listdir("."))


url = "https://www.hitran.org/data/Aerosols/Aerosols-2024/Aerosol_Readme_2024.pdf"
r = requests.get(url)
print(f"Status Code: {r.status_code}")
print(f"Content Snippet: {r.content[:200]}")

print(f"Server: {r.headers.get('Server')}")
print(f"Cache Status: {r.headers.get('X-Cache-Status')}")
