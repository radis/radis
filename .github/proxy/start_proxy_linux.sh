#!/bin/bash

openssl req -x509 -nodes -days 365 -newkey rsa:2048 -subj "/CN=www.hitran.org" -keyout cert.key -out cert.pem -addext "subjectAltName = DNS:www.hitran.org,DNS:hitran.org,DNS:localhost,IP:127.0.0.1"

export NGINX_CACHE_DIR=/home/runner/work/radis/radis/nginx_cache
envsubst '$NGINX_CACHE_DIR' < .github/proxy/nginx.conf > .github/proxy/nginx.conf.gen
sudo nginx -c ${{ github.workspace }}/.github/proxy/nginx.conf.gen -g "user $(whoami);"&

echo "127.0.0.1 www.hitran.org  hitran.org" | sudo tee -a /etc/hosts

sudo cp "${{ github.workspace }}/cert.pem" /usr/local/share/ca-certificates/nginx-proxy.crt
sudo update-ca-certificates
echo "REQUESTS_CA_BUNDLE=/etc/ssl/certs/ca-certificates.crt" >> $GITHUB_ENV
