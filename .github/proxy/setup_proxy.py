# 6. Finally, reroute the www.hitran.org domain to the localhost so requests are sent to the proxy
import platform
import sys

current_os = platform.system()

if current_os == 'Linux':
    HOSTS_PATH = '/etc/hosts'
elif current_os == 'Windows':
    HOSTS_PATH = 'C:\Windows\System32\drivers\etc\hosts'
elif current_os == 'Darwin':
    HOSTS_PATH = '/etc/hosts'
else:
    print(f"Unsupported operating system: {current_os}")
    sys.exit(1)



def append_to_hosts(ip_address, domain_names):
    
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
        print(f"Successfully added '{new_entry[1:]}'")

    except PermissionError:
        print("Permission Denied! This script must be run with admin/root privileges.")
        print(f"   Current OS detected: {current_os}")
        sys.exit(1)


from datetime import datetime, timedelta, timezone
from cryptography import x509
from cryptography.x509.oid import NameOID
from cryptography.hazmat.primitives import hashes
from cryptography.hazmat.primitives.asymmetric import rsa
from cryptography.hazmat.primitives import serialization
from ipaddress import ip_address

def generate_self_signed_cert():
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

    print("Successfully generated cert.key and cert.pem!")


def append_self_signed_cert(self_signed_cert='cert.pem'):
    import certifi
    
    with open(self_signed_cert, 'r') as fr:
        ss_cert_content = fr.read()
    
    print('Certifi:  ',certifi.where())
    print('Requests: ',requests.certs.where())
    
    certifi_fname = certifi.where()
    
    try:
        with open(certifi_fname, 'a') as fa:
            
            fa.write(ss_cert_content)
        print('appended key to certifi!')
    except PermissionError:
        print('Could not append certificates!')
        #TODO: if this happens, make a local copy of the combined certificates and set environment variables REQUESTS_CA_BUNDLE and SSL_CERT_FILE to point at this new file.
        sys.exit(1)
        

if __name__ == "__main__":
    generate_self_signed_cert()
    append_self_signed_cert()
    append_to_hosts("127.0.0.1", ['www.hitran.org', 'hitran.org'])
