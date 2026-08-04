# 6. Finally, reroute the www.hitran.org domain to the localhost so requests are sent to the proxy
import platform
import sys

current_os = platform.system()

if current_os == 'Linux':
    HOSTS_PATH = r'/etc/hosts'
    NGINX_FOLDERS = []
elif current_os == 'Windows':
    HOSTS_PATH = r'C:\Windows\System32\drivers\etc\hosts'
    NGINX_FOLDERS = ['logs', 
                     'temp', 
                     'temp/client_body_temp',
                     'temp/client_body_temperature',
                     'temp/fastcgi_temp',
                     'temp/proxy_remp',
                     'temp/scgi_temp',
                     'temp/uwsgi_temp',
                     ]
elif current_os == 'Darwin':
    HOSTS_PATH = r'/etc/hosts'
    NGINX_FOLDERS = []
else:
    print(f"Unsupported operating system: {current_os}")
    sys.exit(1)



def append_to_hosts(ip_address, domain_names):
    print('-> Adding ', domain_names, 'to hosts... ', end='')
    new_entry = f"\n{ip_address}"
    
    if isinstance(domain_names, str):
        new_entry += f" {domain_names}"
    else:
        for domain_name in domain_names:
            new_entry += f" {domain_name}"
    
    try:
        # 3. Open the file in append mode and write the entry
        with open(HOSTS_PATH, 'a') as hosts_file:
            hosts_file.write(new_entry)
        print('Done!')

    except PermissionError:
        print(f"FAIL - Permission Denied!\nThis script must be run with admin/root privileges (OS={current_os}).")
        return 1


from datetime import datetime, timedelta, timezone
from cryptography import x509
from cryptography.x509.oid import NameOID
from cryptography.hazmat.primitives import hashes
from cryptography.hazmat.primitives.asymmetric import rsa
from cryptography.hazmat.primitives import serialization
from ipaddress import ip_address

def generate_self_signed_cert():
    print('-> Generating self singed certificate... ', end='')
    private_key = rsa.generate_private_key(
        public_exponent=65537,
        key_size=2048
    )

    subject = issuer = x509.Name([
        x509.NameAttribute(NameOID.COMMON_NAME, "www.hitran.org"),
    ])

    cert_builder = (
        x509.CertificateBuilder()
        .subject_name(subject)
        .issuer_name(issuer)
        .public_key(private_key.public_key())
        # Generate a unique serial number
        .serial_number(x509.random_serial_number())
        # Validity timestamps (-days 365)
        .not_valid_before(datetime.now(timezone.utc))
        .not_valid_after(datetime.now(timezone.utc) + timedelta(days=365))
    )

    san_extension = x509.SubjectAlternativeName([
        x509.DNSName("www.hitran.org"),
        x509.DNSName("hitran.org"),
        x509.DNSName("localhost"),
        x509.IPAddress(ip_address("127.0.0.1")),
    ])
    
    cert_builder = cert_builder.add_extension(san_extension, critical=False)
    cert = cert_builder.sign(private_key, hashes.SHA256())

    with open("cert.key", "wb") as f:
        f.write(
            private_key.private_bytes(
                encoding=serialization.Encoding.PEM,
                format=serialization.PrivateFormat.TraditionalOpenSSL,
                encryption_algorithm=serialization.NoEncryption() 
            )
        )

    with open("cert.pem", "wb") as f:
        f.write(cert.public_bytes(serialization.Encoding.PEM))

    print("Done!")


def append_self_signed_cert(cert_file='cert.pem'):
    import certifi
    import requests
    import os
    
    #print('Certifi:  ',certifi.where())
    #print('Requests: ',requests.certs.where())
    
    certs_path = certifi.where()
    
    print(f'-> Appending self-signed certificate to {certs_path}... ', end='') 

    
    # Create a new combined file in the CWD
    with open('ca-certificates.crt', 'w') as fw:
        with open(certs_path, 'r') as fr:
                fw.write(fr.read())
                
        fw.write('\n')
                
        with open(cert_file, 'r') as fr:
            fw.write(fr.read())
   
    print('Done!')
   

def generate_nginx_config(conf_path='.github/proxy/nginx_template.conf'):
    import os
    
    print('-> Generating NGINX config... ', end='')
    
    with open(conf_path, 'r') as fr:
        conf_file = fr.read()
    
    base_dir = os.path.abspath(os.getcwd())
    cache_dir = os.path.join(base_dir, 'nginx_cache')
    
    new_conf_file = conf_file.replace('${NGINX_CACHE_DIR}', cache_dir)

    with open('nginx.conf','w') as fw:
        fw.write(new_conf_file)
        
    for folder in NGINX_FOLDERS:
        os.makedirs(folder, exist_ok=True)
        
    print('Done!')

if __name__ == "__main__":
    generate_self_signed_cert()
    append_self_signed_cert()
    append_to_hosts("127.0.0.1", ['www.hitran.org', 'hitran.org'])
    generate_nginx_config()
