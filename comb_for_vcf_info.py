# Find vcf in other folders based on variant ID from provided CSV file

import os
import csv
import sys

def find_vcf_files(base_dir, target_names):
    """Recursively find VCF files matching target_names (set of IDs)."""
    matches = {}
    for root, _, files in os.walk(base_dir):
        for f in files:
            if f.endswith(".vcf"):
                name = os.path.splitext(f)[0]
                if name in target_names:
                    matches[name] = os.path.join(root, f)
    return matches


def extract_vcf_data(vcf_path):
    """Extract CHROM, POS, ID, REF, ALT from a VCF file."""
    results = []
    with open(vcf_path, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) >= 5:
                chrom, pos, vid, ref, alt = parts[:5]
                results.append([chrom, pos, vid, ref, alt])
    return results


def main(input_csv, base_dir, output_csv):
    # Read variant IDs from CSV (assumes one per row or in first column)
    variant_ids = set()
    with open(input_csv, "r") as f:
        reader = csv.reader(f)
        for row in reader:
            if row:
                variant_ids.add(row[0].strip())

    print(f"Loaded {len(variant_ids)} variant IDs")

    # Find matching VCF files
    vcf_files = find_vcf_files(base_dir, variant_ids)
    print(f"Found {len(vcf_files)} matching VCF files")

    # Extract and write output
    with open(output_csv, "w", newline="") as out:
        writer = csv.writer(out)
        writer.writerow(["chrom", "pos", "id", "ref", "alt"])

        for vid, vcf_path in vcf_files.items():
            records = extract_vcf_data(vcf_path)
            for r in records:
                writer.writerow(r)

    print(f"Output written to {output_csv}")


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script.py <input_csv> <vcf_directory> <output_csv>")
        sys.exit(1)

    input_csv = sys.argv[1]
    base_dir = sys.argv[2]
    output_csv = sys.argv[3]

    main(input_csv, base_dir, output_csv)