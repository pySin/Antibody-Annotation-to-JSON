# Description: a script that parses and converts to JSON an annotation format (designed by Professor Andrew Martin) for antibody-based drugs


# import sys
import os
import re
import json
import argparse

# Regular expression patterns to extract specific data
p_request = re.compile(r'^Request:\s*(\d+);', re.MULTILINE)
p_antigen = re.compile(r'Antigen:\s*(.+?)\s*\((Gene[a-zA-Z0-9]*)\);', re.MULTILINE)
p_note = re.compile(r'Note:\s*(.+?);', re.MULTILINE)
p_chain = re.compile(r'(Heavy|Light) Chain\[(\d+)\]:\s*([A-Z]+);', re.MULTILINE)

# Function to parse a single file and convert it to JSON
def txt_to_json(txtpath, outpath):
    # Dictionary to store annotations
    annotations = {}

    # Open and read the file
    with open(txtpath, 'r') as f:
        txt = f.read()

    # Extract the 'Request' field
    request_match = p_request.search(txt)
    if request_match:
        annotations['Request'] = request_match.group(1)
    else:
        annotations['Request'] = "Unknown"

    # Extract the 'Antigen' field
    antigen_match = p_antigen.search(txt)
    if antigen_match:
        annotations['Antigen'] = {
            "Name": antigen_match.group(1),
            "Gene": antigen_match.group(2)
        }

    # Extract all 'Notes'
    notes = p_note.findall(txt)
    if notes:
        annotations['Notes'] = notes

    # Extract 'Chains' (e.g., Heavy and Light chains)
    chain_matches = p_chain.findall(txt)
    if chain_matches:
        annotations['Chains'] = []
        for chain_type, instance, sequence in chain_matches:
            annotations['Chains'].append({
                'Type': chain_type,
                'Instance': instance,
                'Sequence': sequence
            })

    # Write the structured data to a JSON file
    with open(outpath, 'w') as out:
        json.dump(annotations, out, indent=4)

    print(f"Processed file: {txtpath} -> {outpath}")

# Process a single file
def parse_file(filepath, outdir):
    outpath = os.path.join(outdir, os.path.basename(filepath).replace(".txt", ".json"))
    txt_to_json(filepath, outpath)

# Process all files in a directory
def parse_directory(directory, outdir):
    for filename in os.listdir(directory):
        if filename.endswith(".txt"):
            filepath = os.path.join(directory, filename)
            parse_file(filepath, outdir)

# Main function to handle command-line arguments
def main():
    parser = argparse.ArgumentParser(description="Parse antibody annotation files to JSON.")
    parser.add_argument(
        "-f", "--file",
        type=str,
        help="Path to a single text file to process."
    )
    parser.add_argument(
        "-d", "--directory",
        type=str,
        help="Path to a directory containing text files to process."
    )
    parser.add_argument(
        "-o", "--output",
        type=str,
        default="json",
        help="Output directory for JSON files (default: ./json)."
    )
    args = parser.parse_args()

    # Create the output directory if it does not exist
    if not os.path.exists(args.output):
        os.makedirs(args.output)

    # Handle input data
    if args.file:
        parse_file(args.file, args.output)
    elif args.directory:
        parse_directory(args.directory, args.output)
    else:
        print("Error: Either a file or directory path must be specified.")
        parser.print_help()

if __name__ == "__main__":
    main()

# Usage:
# - To process a single file:
#   python parser.py -f path/to/file.txt -o path/to/output_directory
#
# - To process all files in a directory:
#   python parser.py -d path/to/directory -o path/to/output_directory
#
# - If no output directory is specified, JSON files will be saved in the default 'json/' directory.
#
# Examples:
# - python parser.py -f samples/example.txt
# - python parser.py -d samples/ -o output_json


# (DONE) TODO: flags/options for different types of user input: file path (one or many) vs directory path
# TODO: split some of the records into two e.g. Antigen-Gene for the gene name in parentheses, and CDRKabatH1-Range for the residue range in parentheses
