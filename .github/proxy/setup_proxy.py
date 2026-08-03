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
        print(f"Successfully added '{new_entry}'")

    except PermissionError:
        print("Permission Denied! This script must be run with admin/root privileges.")
        print(f"   Current OS detected: {current_os}")
        sys.exit(1)


if __name__ == "__main__":
    append_to_hosts("127.0.0.1", ['www.hitran.org', 'hitran.org'])
