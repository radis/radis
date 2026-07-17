# import os
import requests
import urllib3
# import sys 

# print("CA_bundle:", os.environ["REQUESTS_CA_BUNDLE"])
# print("${{ github.workspace }} path:", os.path.abspath("."))
# print("${{ github.workspace }} contents:", os.listdir("."))

# cert_path = sys.argv[-1]
# print('CERT_FILE', cert_path)
# # 1. Direct OpenSSL and requests to trust your CA bundle
# os.environ["SSL_CERT_FILE"] = cert_path
# os.environ["REQUESTS_CA_BUNDLE"] = cert_path

# # 2. Tell Python where your Nginx proxy is located
# os.environ["HTTP_PROXY"] = "http://127.0.0.1:443"
# os.environ["HTTPS_PROXY"] = "http://127.0.0.1:443"

# # 3. Handle the routing logic
# # 'requests' reads the wildcard '*' perfectly.
# # For 'urllib', explicitly list common domains to bypass instead (e.g., '.com', '.org')
# # os.environ["no_proxy"] = "*, localhost, 127.0.0.1, .com, .org, .fr"
# os.environ["NO_PROXY"] = "*, localhost, 127.0.0.1, .com, .org, .fr"


headers = {
    'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36'
}

url1 = "https://www.hitran.org/data/Aerosols/Aerosols-2024/Aerosol_Readme_2024.pdf"
url2 = "https://hitran.org/data/Aerosols/Aerosols-2024/Aerosol_Readme_2024.pdf"
url3 = "https://geisa.aeris-data.fr/wp-content-aeris/uploads/sites/33/2021/07/table_abundance_GEISA_2020.pdf"

for i, url in enumerate([url1, url2, url3]):

        
    print(i+1, 'requests')
    r = requests.get(url)
    print(f"Status Code: {r.status_code}")
    print(f"Content Snippet: {r.content[:200]}")
    print(f"Server: {r.headers.get('Server')}")
    print(f"Cache Status: {r.headers.get('X-Cache-Status')}")
    

    print('\n')
    
    print(i+1, 'urllib3')
    http = urllib3.PoolManager()
    r = http.request("GET", url, headers=headers)
    print(f"Status Code: {r.status}")
    print(f"Content Snippet: {r.data[:200]}")
    print(f"Server: {r.headers.get('Server')}")
    print(f"Cache Status: {r.headers.get('X-Cache-Status')}")

    

    print('\n=====================\n')
