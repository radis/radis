openssl req -x509 -nodes -days 365 -newkey rsa:2048 -subj "/CN=www.hitran.org" -keyout cert.key -out cert.pem -addext "subjectAltName = DNS:www.hitran.org,DNS:hitran.org,DNS:localhost,IP:127.0.0.1"

cd C:\tools\nginx-*\
Stop-Service -Name "W3SVC" -Force
$env:NGINX_CACHE_DIR = "D:/a/radis/radis/nginx_cache"
(Get-Content ${{ github.workspace }}/.github/proxy/nginx.conf) -replace '\$\{NGINX_CACHE_DIR\}', $env:NGINX_CACHE_DIR | Set-Content ${{ github.workspace }}/.github/proxy/nginx.conf.gen
Start-Process ".\nginx.exe" -ArgumentList "-c", "${{ github.workspace }}/.github/proxy/nginx.conf.gen"
cd ${{ github.workspace }}

Add-Content -Path C:\Windows\System32\drivers\etc\hosts -Value "`r`n127.0.0.1 www.hitran.org hitran.org"

Import-Certificate -FilePath .\cert.pem -CertStoreLocation Cert:\LocalMachine\Root
Copy-Item -Path .\cert.pem -Destination "${{ github.workspace }}/combined_cacert.pem"

pip install certifi
$pythonCertPath = python -c "import certifi; print(certifi.where())"
Get-Content -Path $pythonCertPath | Add-Content -Path "${{ github.workspace }}/combined_cacert.pem"

"REQUESTS_CA_BUNDLE=${{ github.workspace }}/combined_cacert.pem" | Out-File -FilePath $env:GITHUB_ENV -Encoding utf8 -Append

Get-Process nginx
