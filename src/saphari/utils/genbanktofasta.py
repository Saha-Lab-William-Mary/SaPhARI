#!/usr/bin/env python

import argparse
import os
from Bio import SeqIO

def process_genbank(input_filename):
    results = []

    for record in SeqIO.parse(input_filename, "genbank"):
        for feature in record.features:
            if feature.type == "CDS":
                try:
                    protein_name = feature.qualifiers.get("product", [""])[0].replace(" ", "_")  # Replace spaces with underscores
                    protein_id = feature.qualifiers.get('protein_id', [""])[0]
                    start = feature.location.start
                    end = feature.location.end
                    translation = feature.qualifiers.get('translation', [None])[0]

                    # Skip if translation is missing
                    if translation is None:
                        continue

                    # Prepare the header
                    header = f">{protein_name}|{start}|{end}"
                    if protein_id:  # Add protein ID only if it's present
                        header += f"|{protein_id}"

                    fasta_entry = f"{header}\n{translation}"
                    results.append(fasta_entry)
                except KeyError:
                    continue  # Skip if any required information is missing
    
    return results

def main():
    # Setup argument parser
    argparser = argparse.ArgumentParser(description='Process a GenBank file to extract CDS regions and format them as FASTA.')
    argparser.add_argument('--input', type=str, required=True, help='Input GenBank file to process')
    ARGS = argparser.parse_args()

    input_filename = ARGS.input
    output_filename = f"{os.path.splitext(input_filename)[0]}_gbtofasta_{os.path.splitext(input_filename)[1]}"

    # Process the input file
    results = process_genbank(input_filename)

    # Write the results to the output file
    with open(output_filename, 'w') as outfile:
        for entry in results:
            outfile.write(f"{entry}\n")

if __name__ == '__main__':
    main()
