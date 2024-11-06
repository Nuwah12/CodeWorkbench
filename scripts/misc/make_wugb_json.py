#!/usr/bin/env python3

"""
This script scans a web server's directory (based on `servers` dict),
and saves a JSON file that includes all supported data files in this
directory.

The output JSON file can be loaded directly into Washington University
Genome Browser (WUGB).

Two arguments are required by the script:
  * genome type, such as `hg19` or `hg38`
  * directory name

Two new files will be generated in the input directory:
  * `wugb.json`: JSON file for WUGB
  * `wugb_url.txt`: the URL of WUGB (also shown at the end of this program)
"""

import json
import os
import socket
import sys

# Web servers config
servers = {
    # plutus
    'plutus': {
        'www_root': '/mnt/data1/www/html/',
        'main_url': 'https://plutus.faryabilab.com/',
    },

    # simurgh
    'simurgh': {
        'www_root': '/mnt/data0/www/html/',
        'main_url': 'https://simurgh.faryabilab.com/',
    },
}

# Types of data that will be rendered in WashU Genome Browser
data_types = {
    '.bw': 'bigwig',
    '.cool': 'cool',
    '.hic': 'hic',
}

# Main URL of WUGB
WUGB_URL = "https://epigenomegateway.wustl.edu/browser/"

# Output filenames
JSON_FILENAME = "wugb.json"
URL_FILENAME = "wugb_url.txt"

def syntax():
    print("Syntax:")
    print("  make_wugb_json.py [genome_type] [data_directory]\n")

    print("For example:")
    print("  make_wugb_json.py hg38 my_data")


def chk_data_dir(dir_name, data_root):
    """
    Check whether:
     * `dir_name` is a directory
     * `dirname` is a sub-directory of `data_root`

    Return the absolute path of `dir_name`.
    """

    abs_path = os.path.abspath(dir_name)
    if not os.path.isdir(abs_path):
        print(f"ERROR: '{dir_name}' is not a directory")
        sys.exit(2)

    if not abs_path.startswith(data_root):
        print(f"ERROR: '{dir_name}' is not located in '{data_root}'")
        sys.exit(3)

    return abs_path


# Main
if __name__ == '__main__':
    if len(sys.argv) != 3:
        syntax()
        sys.exit(1)

    # Exit if the server is not found in `servers`
    server_name = socket.gethostname()
    if server_name not in servers:
        print(f"ERROR: '{server_name}' not supported")
        sys.exit(3)

    www_root = servers[server_name]['www_root']
    main_url = servers[server_name]['main_url']

    genome_type = sys.argv[1]
    data_dir = sys.argv[2]

    abs_data_dir = chk_data_dir(data_dir, www_root)
    wugb_hub = list()
    for root, _, files in os.walk(abs_data_dir):
        # Make sure `root` is accessible by all users
        os.chmod(root, 0o755)

        for f in files:
            # Skip the files whose types are not supported
            _, f_type = os.path.splitext(f)
            if f_type not in data_types:
                continue

            abs_file_path = os.path.join(root, f)
            sub_url = os.path.relpath(abs_file_path, start=www_root)
            data_url = main_url + sub_url
            hub_entry = dict()
            hub_entry["type"] = data_types[f_type]
            hub_entry["url"] = data_url
            hub_entry["name"] = f.split('.')[0]
            hub_entry["showOnHubLoad"] = True

            # Special options for `cool` and `hic` data types, see:
            # https://epigenomegateway.readthedocs.io/en/latest/datahub.html#example-hic-track
            if f_type in ['cool', 'hic']:
                hub_entry['options'] = {
                    'displayMode': 'arc',
                }

            wugb_hub.append(hub_entry) # add entry to hub

    # Create `wugb.json` in `data_dir`
    json_path = os.path.join(abs_data_dir, JSON_FILENAME)
    with open(json_path, "w") as ofh:
        json.dump(wugb_hub, ofh, indent=2)
        ofh.write("\n")  # end the JSON file with a newline character

    json_url = main_url + os.path.relpath(json_path, start=www_root)

    # Create `wugb_url.txt` in `bw_dir`:
    url_path = os.path.join(abs_data_dir, URL_FILENAME)
    url_str = f"{WUGB_URL}?genome={genome_type}&hub={json_url}"
    with open(url_path, "w") as ofh:
        ofh.write(f"{url_str}\n")

    # Print `json_path` and `url_str` on stdout
    print(f"JSON file is saved as: {json_path}")
    print(f"View it at: {url_str}")
    print(f"URL is saved as: {url_path}")
