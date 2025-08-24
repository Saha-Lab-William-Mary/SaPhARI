#!/usr/bin/env python

import argparse
import os
from collections import Counter
import re

def process_file(input_filename):
    orf_to_descriptions = {}
    skipped_lines = []

    # Read data from the file
    with open(input_filename, 'r') as infile:
        lines = infile.readlines()
        for line in lines:
            parts = line.strip().split('\t')
            if len(parts) < 2:
                skipped_lines.append(line.strip())
                continue  # Skip lines that don't have at least 2 columns
            
            # Extract ORF from the part before the first '|'
            orf = parts[0].split('|')[0]
            
            # Extract the description part
            description_match = re.search(r'([^\[]+)', parts[1])
            if description_match:
                description = description_match.group(1).strip()
                if orf not in orf_to_descriptions:
                    orf_to_descriptions[orf] = []
                orf_to_descriptions[orf].append((description, line.strip()))
            else:
                skipped_lines.append(line.strip())

    orf_to_common_description = {}
    
    # Process each ORF group
    for orf, descriptions in orf_to_descriptions.items():
        # Separate hypothetical and non-hypothetical descriptions
        non_hypothetical = [desc for desc in descriptions if 'hypothetical' not in desc[0].lower()]
        hypothetical = [desc for desc in descriptions if 'hypothetical' in desc[0].lower()]

        if non_hypothetical:
            # Find the most common description among non-hypothetical descriptions
            descriptions_text = [desc[0] for desc in non_hypothetical]
            most_common_description = Counter(descriptions_text).most_common(1)[0][0]
            # Find the first line containing the most common description
            for desc, line in descriptions:
                if desc == most_common_description:
                    orf_to_common_description[orf] = line
                    break
        else:
            # If all lines are hypothetical, just take the first hypothetical line
            if hypothetical:
                orf_to_common_description[orf] = hypothetical[0][1]

    # Handle any missing ORFs that might have been skipped
    for orf, descriptions in orf_to_descriptions.items():
        if orf not in orf_to_common_description:
            orf_to_common_description[orf] = descriptions[0][1]

    # Gather the results
    results = list(orf_to_common_description.values())
    results.extend(skipped_lines)
    
    return results

def main():
    # Setup argument parser
    argparser = argparse.ArgumentParser()
    argparser.add_argument('--input', type=str, required=True, help='Input file to process')
    ARGS = argparser.parse_args()

    input_filename = ARGS.input
    output_filename = f"{os.path.splitext(input_filename)[0]}_filtered{os.path.splitext(input_filename)[1]}"

    # Process the input file
    results = process_file(input_filename)

    # Write the results to the output file
    with open(output_filename, 'w') as outfile:
        for line in results:
            outfile.write(f"{line}\n")

if __name__ == '__main__':
    main()
