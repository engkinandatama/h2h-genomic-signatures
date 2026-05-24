#!/usr/bin/env python3
"""
Diagnostic script: Dump all unique host and location values from NCBI ZIP files.
Usage: python inspect_metadata.py <path_to_dataset.zip>
"""
import sys
import json
import zipfile
from collections import Counter

def inspect_zip(zip_path):
    print(f"\n{'='*60}")
    print(f"Inspecting: {zip_path}")
    print(f"{'='*60}")

    with zipfile.ZipFile(zip_path, 'r') as z:
        metadata_files = [f for f in z.namelist() if f.endswith('data_report.jsonl')]
        if not metadata_files:
            print("ERROR: data_report.jsonl not found!")
            return

        host_counter = Counter()
        geo_counter = Counter()
        null_host_count = 0
        total = 0

        with z.open(metadata_files[0]) as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                try:
                    record = json.loads(line.decode('utf-8'))
                except:
                    continue

                total += 1

                # Extract host
                host_obj = record.get('host')
                if not isinstance(host_obj, dict):
                    null_host_count += 1
                    host_counter['<NO HOST FIELD>'] += 1
                else:
                    host_name = host_obj.get('organismName') or host_obj.get('name') or '<EMPTY>'
                    host_counter[host_name] += 1

                # Extract location
                loc_obj = record.get('location')
                if not isinstance(loc_obj, dict):
                    geo_counter['<NO LOCATION FIELD>'] += 1
                else:
                    geo = loc_obj.get('geographicLocation') or '<EMPTY>'
                    geo_counter[geo] += 1

        print(f"\nTotal records: {total}")
        print(f"Records with NO host field: {null_host_count}")

        print(f"\n--- TOP 30 HOST VALUES ---")
        for host, count in host_counter.most_common(30):
            print(f"  {count:5d}x  '{host}'")

        print(f"\n--- TOP 30 GEOGRAPHIC LOCATIONS ---")
        for geo, count in geo_counter.most_common(30):
            print(f"  {count:5d}x  '{geo}'")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python inspect_metadata.py <dataset.zip> [dataset2.zip ...]")
        sys.exit(1)
    for path in sys.argv[1:]:
        inspect_zip(path)
