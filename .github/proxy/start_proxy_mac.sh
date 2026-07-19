#!/bin/bash

openssl req -x509 -nodes -days 365 -newkey rsa:2048 -subj "/CN=www.hitran.org" -keyout cert.key -out cert.pem -addext "subjectAltName = DNS:www.hitran.org,DNS:hitran.org,DNS:localhost,IP:127.0.0.1"

brew install nginx
export NGINX_CACHE_DIR=/Users/runner/work/radis/radis/nginx_cache
envsubst '$NGINX_CACHE_DIR' < .github/proxy/nginx.conf > .github/proxy/nginx.conf.gen
sudo nginx -c ${{ github.workspace }}/.github/proxy/nginx.conf.gen -g "user $(whoami) staff;"&

echo "127.0.0.1 www.hitran.org hitran.org" | sudo tee -a /etc/hosts

sudo security add-trusted-cert -d -r trustRoot -k /Library/Keychains/System.keychain cert.pem

python3 -c "import certifi; print(certifi.where())" > certifi_path.txt
cp $(cat certifi_path.txt) combined_ca.pem
cat cert.pem >> combined_ca.pem

echo "REQUESTS_CA_BUNDLE=${{ github.workspace }}/combined_ca.pem" >> $GITHUB_ENV
echo "SSL_CERT_FILE=${{ github.workspace }}/combined_ca.pem" >> $GITHUB_ENV
