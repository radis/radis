# 6. Finally, reroute the www.hitran.org domain to the localhost so requests are sent to the proxy
import platform
import sys


def append_to_hosts(ip_address, domain_name):
    # 1. Determine the path to the hosts file based on the OS
    current_os = platform.system()

    if current_os == "Windows":
        hosts_path = r"C:\Windows\System32\drivers\etc\hosts"
    elif current_os in ["Linux", "Darwin"]:  # Darwin is macOS
        hosts_path = "/etc/hosts"
    else:
        print(f"Unsupported operating system: {current_os}")
        sys.exit(1)

    # 2. Format the line to append
    new_entry = f"\n{ip_address} {domain_name}\n"

    try:
        # 3. Open the file in append mode and write the entry
        with open(hosts_path, "a") as hosts_file:
            hosts_file.write(new_entry)
        print(f"Successfully added '{domain_name}' to {hosts_path}")

    except PermissionError:
        print("Permission Denied! This script must be run with admin/root privileges.")
        print(f"   Current OS detected: {current_os}")
        sys.exit(1)


if __name__ == "__main__":
    append_to_hosts("127.0.0.1", "www.hitran.org")
