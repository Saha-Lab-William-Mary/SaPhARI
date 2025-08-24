#!/usr/bin/env python3
import sys
import os
import re
import pandas as pd

def parse_gene_line(line):
    """
    Parse a gene line of the form:
      <gene_id>|<start>|<stop>|<strand>|[optional extra fields...]
    Returns (gene_id, start, stop, strand) or None if malformed.
    """
    parts = line.strip().split('|')
    if len(parts) < 4:
        return None
    gene_id = parts[0].strip()
    start = parts[1].strip()
    stop = parts[2].strip()
    strand = parts[3].strip()
    return gene_id, start, stop, strand

def parse_blast_line(line):
    """
    Parse a BLAST‐style line. Returns (protein_id, protein_name,
    percent_identity, percent_coverage, blast_evalue).
    """
    line = line.strip()
    tokens = re.split(r'\t+', line)
    if len(tokens) >= 4:
        first_field = tokens[0].strip()
        parts = first_field.split(" ", 1)
        protein_id = parts[0]
        protein_name = parts[1].strip() if len(parts) > 1 else ""
        percent_identity = tokens[1].strip()
        percent_coverage = tokens[2].strip()
        blast_evalue = tokens[3].strip()
    else:
        tokens = line.split()
        if len(tokens) < 5:
            return "", "", "", "", ""
        protein_id = tokens[0]
        percent_identity = tokens[-3]
        percent_coverage = tokens[-2]
        blast_evalue = tokens[-1]
        protein_name = " ".join(tokens[1:-3])
    return protein_id, protein_name, percent_identity, percent_coverage, blast_evalue

def process_file(filename):
    """
    Read a Satellite‐output “*_output.txt” file and return a list of region dicts:
      [
        {
          "name": "PICI REGION 1 - file_basename",
          "upstream": [<gene_block>, ...],
          "core":     [<gene_block>, ...],
          "downstream":[<gene_block>, ...]
        },
        ...
      ]
    Each gene_block is a dict with keys:
      Region, Section, ORF No., Start, Stop, Strand,
      Protein ID, %ID, %Cov, Blast e-value, Protein Name
    """
    regions = []
    current_region = None
    current_section = None
    region_has_core = False
    current_file_basename = None

    with open(filename, 'r') as f:
        lines = f.readlines()

    i = 0
    while i < len(lines):
        line = lines[i].strip()

        # Skip blank lines
        if not line:
            i += 1
            continue

        # Capture "Matching lines for family 'FAM' in file '...filename...'":
        m = re.match(r"Matching lines for family '.+' in file '(.+)'", line)
        if m:
            full_path = m.group(1)
            fname = os.path.basename(full_path)
            fname_noext = os.path.splitext(fname)[0]
            # Strip any "_prot_orfed" suffix if present
            if "_prot_orfed" in fname_noext:
                current_file_basename = fname_noext.split("_prot_orfed")[0]
            else:
                current_file_basename = fname_noext
            i += 1
            continue

        # Detect "<FAMILY> REGION n:" header
        if re.search(r"REGION\s+\d+:$", line):
            if current_region is not None:
                regions.append(current_region)
            region_header = line.rstrip(":")  # e.g. "PICI REGION 1"
            # Combine region header with file basename
            if current_file_basename:
                region_name = f"{region_header} - {current_file_basename}"
            else:
                region_name = region_header
            current_region = {
                "name": region_name,
                "upstream": [],
                "core": [],
                "downstream": []
            }
            current_section = None
            region_has_core = False
            i += 1
            continue

        # Skip lines of "===="
        if re.fullmatch(r"[=]+", line):
            i += 1
            continue

        # Detect upstream flanking vs. downstream flanking headers
        if line.startswith("Following Flanking Genes:"):
            current_section = "Downstream"
            i += 1
            continue
        if line.startswith("Flanking Genes:"):
            current_section = "Upstream" if not region_has_core else "Downstream"
            i += 1
            continue

        # "Core Protein Region:" marks the Core block
        if line.startswith("Core Protein Region:"):
            current_section = "Core"
            region_has_core = True
            i += 1
            continue

        # If we're inside a section and this line has a '|' (gene line), parse it
        if current_section and "|" in line:
            gene_line = line
            # Next non‐empty, non‐header line is the BLAST line
            j = i + 1
            blast_line = ""
            while j < len(lines):
                next_line = lines[j].strip()
                if not next_line:
                    j += 1
                    continue
                if next_line.endswith("Genes:") or next_line.endswith("Region:"):
                    j += 1
                    continue
                blast_line = next_line
                break

            parsed_gene = parse_gene_line(gene_line)
            if not parsed_gene:
                i = j + 1
                continue

            gene_id, start, stop, strand = parsed_gene
            orf_num = gene_id.split("_")[-1] if "_" in gene_id else gene_id
            prot_id, prot_name, perc_id, perc_cov, blast_evalue = parse_blast_line(blast_line)

            gene_block = {
                "Region": current_region["name"],
                "Section": current_section,
                "ORF No.": orf_num,
                "Start": start,
                "Stop": stop,
                "Strand": strand,
                "Protein ID": prot_id,
                "%ID": perc_id,
                "%Cov": perc_cov,
                "Blast e-value": blast_evalue,
                "Protein Name": prot_name
            }

            if current_section == "Core":
                current_region["core"].append(gene_block)
            elif current_section == "Upstream":
                current_region["upstream"].append(gene_block)
            elif current_section == "Downstream":
                current_region["downstream"].append(gene_block)

            i = j + 1
            continue

        # Otherwise, skip unrecognized lines
        i += 1

    # Append the last region if present
    if current_region is not None:
        regions.append(current_region)
    return regions

def export_regions_to_tsv(regions, output_tsv):
    """
    Given a list of region dicts (as returned by process_file), flatten them and
    write to a TSV with columns:
      Region, Section, ORF No., Start, Stop, Strand,
      Protein ID, %ID, %Cov, Blast e-value, Protein Name
    """
    all_blocks = []
    for region in regions:
        all_blocks.extend(region["upstream"])
        all_blocks.extend(region["core"])
        all_blocks.extend(region["downstream"])

    cols = [
        "Region", "Section", "ORF No.", "Start", "Stop", "Strand",
        "Protein ID", "%ID", "%Cov", "Blast e-value", "Protein Name"
    ]

    if not all_blocks:
        df_empty = pd.DataFrame(columns=cols)
        df_empty.to_csv(output_tsv, sep="\t", index=False)
        return

    df = pd.DataFrame(all_blocks)
    df = df[cols]
    df.to_csv(output_tsv, sep="\t", index=False)

def main():
    if len(sys.argv) < 2:
        print("Usage: python convert_output.py input_file.txt [output_file.tsv]")
        sys.exit(1)

    input_file = sys.argv[1]
    if len(sys.argv) >= 3:
        output_tsv = sys.argv[2]
    else:
        base = os.path.splitext(os.path.basename(input_file))[0]
        output_tsv = base + ".tsv"

    regions = process_file(input_file)
    export_regions_to_tsv(regions, output_tsv)
    print(f"Exported TSV to {output_tsv}")

if __name__ == "__main__":
    main()
