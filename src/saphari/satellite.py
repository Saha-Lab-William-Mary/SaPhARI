from Bio import SeqIO
import os

class Satellite:
    
    def __init__(self):
        self.families = []

    def add_family(self, title, proteins=[], length=None, minnumber=None, forbidden=[], specific_distances=None):
        new_family = self.Family(title, proteins, length, minnumber, forbidden, specific_distances)
        self.families.append(new_family)
        return new_family
    
    def get_family(self, title):
        for family in self.families:
            if family.title == title:
                return family
        return None 
    
    def process_files(self, file_list, use_cds, output_type, familiestosearchfor, inputpath, outDir):
        for full_path in file_list:
            if os.path.isfile(full_path):  # Ensure it's a file and not a directory
                file_size = os.path.getsize(full_path)
                print(f"Checking file: {full_path}")
                print(f"File size: {file_size} bytes")

                if file_size > 0:  # Check if the file has non-zero size
                    print(f"File '{full_path}' is ready.")
                    
                    # Derive base_name by stripping any annotation‐pipeline suffix beginning at "_prot_orfed"
                    fname = os.path.basename(full_path)
                    fname_noext = os.path.splitext(fname)[0]
                    if "_prot_orfed" in fname_noext:
                        base_name = fname_noext.split("_prot_orfed")[0]
                    else:
                        base_name = fname_noext

                    for family in self.families:
                        if family.title in familiestosearchfor:
                            if use_cds:
                                matching_lines = family.find_it_gb_CDS(full_path)
                            else:
                                matching_lines = family.find_it(full_path)

                            # Correct output directory for family
                            family_output_dir = os.path.join(inputpath, outDir, 'RESULTS', family.title)
                            os.makedirs(family_output_dir, exist_ok=True)

                            # Construct the output file path correctly
                            if output_type == 'pfpf':
                                family_output_path = os.path.join(
                                    family_output_dir,
                                    f"{base_name}_{family.title}_output.txt"
                                )
                                open_mode = 'w'
                            else:  # 'pf'
                                family_output_path = os.path.join(
                                    family_output_dir,
                                    f"{family.title}_output.txt"
                                )
                                open_mode = 'a'

                            with open(family_output_path, open_mode) as family_output_file:
                                family_output_file.write(
                                    f"Matching lines for family '{family.title}' in file '{full_path}':\n"
                                )
                                family_output_file.writelines(line + '\n' for line in matching_lines)
                                family_output_file.write('\n')  # Add a newline for better separation between files
                else:
                    print(f"File '{full_path}' is empty. Skipping.")
            else:
                print(f"{full_path} is not a file.")

    class Family:
        def __init__(self, title, proteins=[], length=None, minnumber=None, forbidden=[], specific_distances=None, check_forbidden_flanks=False):
            self.title = title
            self.proteins = proteins
            self.length = length
            self.number = minnumber
            self.forbidden = forbidden
            self.specific_distances = specific_distances if specific_distances else {}
            # New flag: if True, the flanking regions will be checked for forbidden proteins.
            self.check_forbidden_flanks = check_forbidden_flanks

        def set_proteins(self, proteins):
            self.proteins = proteins

        def set_length(self, length):
            self.length = length

        def set_forbidden(self, forbidden):
            self.forbidden = forbidden

        def find_it(self, filepath):
            found_lines = []
            all_regions = []

            with open(filepath, 'r') as file:
                lines = file.readlines()

            def extract_proteins(line):
                line_lower = line.lower()
                proteins_found = set()
                for protein in self.proteins:
                    if isinstance(protein, tuple):
                        if any(p.lower() in line_lower for p in protein):
                            proteins_found.add(protein[0].lower())
                    else:
                        if protein.lower() in line_lower:
                            proteins_found.add(protein.lower())
                return proteins_found

            def contains_forbidden(line):
                line_lower = line.lower()
                for protein in self.forbidden:
                    if protein.lower() in line_lower:
                        return True
                return False

            # Build index mapping line content → index
            line_to_index = {line.strip(): idx for idx, line in enumerate(lines)}

            # Collect all potential regions
            for i, line in enumerate(lines):
                line_parts = line.strip().split('|')
                if len(line_parts) < 3:
                    continue
                start_position = int(line_parts[1])
                
                group_line = []
                protein_counter = set()
                min_stop_position = None

                for linex in lines[i:]:  # Continue till the end of file
                    linex_parts = linex.strip().split('|')
                    if len(linex_parts) < 3:
                        continue
                    stop_position = int(linex_parts[2])

                    if min_stop_position is None:
                        min_stop_position = stop_position

                    # Check window length
                    if (stop_position - start_position <= self.length):
                        if contains_forbidden(linex):
                            break
                        proteins_found = extract_proteins(linex)
                        if proteins_found:
                            group_line.append(linex.strip())
                            protein_counter.update(proteins_found)
                            min_stop_position = stop_position
                    else:
                        # If the current region is too long, break
                        break

                if len(protein_counter) >= self.number and len(group_line) > 0:
                    all_regions.append((set(group_line), i, len(group_line), start_position, min_stop_position))

            # Sort regions by their number of lines in descending order
            all_regions.sort(key=lambda x: x[2], reverse=True)

            # Ensure regions are unique
            unique_regions = []
            seen_regions = set()

            def is_subset(region, existing_regions):
                region_lines = set(region[0])
                for existing in existing_regions:
                    if region_lines.issubset(existing[0]):
                        return True
                return False

            for region in all_regions:
                if not is_subset(region, unique_regions):
                    region_lines = set(region[0])
                    seen_regions.add(frozenset(region_lines))
                    unique_regions.append(region)

            # Add auxiliary functions
            def get_flanking_lines(region_start_index, region_end_index):
                previous_flanking = []
                following_flanking = []

                for idx in range(max(0, region_start_index - 5), region_start_index):
                    previous_flanking.append(lines[idx].strip())
                
                for idx in range(region_end_index + 1, min(len(lines), region_end_index + 6)):
                    following_flanking.append(lines[idx].strip())
                
                return previous_flanking, following_flanking

            def split_line_content(line):
                parts = line.split('\t', 1)
                if len(parts) > 1:
                    return [parts[0].strip(), parts[1].strip()]
                return [line.strip()]

            def check_specific_distances(region_lines):
                # Store the positions of found proteins
                positions = {protein: None for protein_pair in self.specific_distances for protein in protein_pair}

                # Iterate over each specified pair of proteins and their maximum allowable distance
                for (protein1, protein2), max_distance in self.specific_distances.items():
                    found_protein1 = found_protein2 = False
                    for line in region_lines:
                        # Check for protein1
                        if isinstance(protein1, tuple):
                            if any(p.lower() in line.lower() for p in protein1):
                                parts = line.split('|')
                                positions[protein1] = (int(parts[1]), int(parts[2]))  # Store start and stop positions
                                found_protein1 = True
                        else:
                            if protein1.lower() in line.lower():
                                parts = line.split('|')
                                positions[protein1] = (int(parts[1]), int(parts[2]))  # Store start and stop positions
                                found_protein1 = True

                        # Check for protein2
                        if isinstance(protein2, tuple):
                            if any(p.lower() in line.lower() for p in protein2):
                                parts = line.split('|')
                                positions[protein2] = (int(parts[1]), int(parts[2]))  # Store start and stop positions
                                found_protein2 = True
                        else:
                            if protein2.lower() in line.lower():
                                parts = line.split('|')
                                positions[protein2] = (int(parts[1]), int(parts[2]))
                                found_protein2 = True

                        if found_protein1 and found_protein2:
                            pos1_end = positions[protein1][1]
                            pos2_start = positions[protein2][0]
                            distance = abs(pos2_start - pos1_end)
                            if distance > max_distance:
                                return False
                            # Reset for next pair
                            positions[protein1] = positions[protein2] = None
                            found_protein1 = found_protein2 = False
                            continue

                # If all specified pairs meet their distance constraints, return True
                return True

            found_lines = []
            for idx, (region_lines, _, _, region_start_position, region_end_position) in enumerate(unique_regions):
                if not check_specific_distances(region_lines):
                    continue

                group_name = f"{self.title} REGION {idx + 1}"
                found_lines.append(f"\n{group_name}:")
                found_lines.append("=" * len(group_name))

                # Determine indices in the original lines
                region_start_index = min(line_to_index[line] for line in region_lines if line in line_to_index)
                region_end_index = max(line_to_index[line] for line in region_lines if line in line_to_index)

                # Add flanking lines
                previous_flanking, following_flanking = get_flanking_lines(region_start_index, region_end_index)

                if self.check_forbidden_flanks:
                    if any(contains_forbidden(line) for line in previous_flanking + following_flanking):
                        continue

                if previous_flanking:
                    found_lines.append("\nFlanking Genes:")
                    found_lines.append("")
                    for line in previous_flanking:
                        split_lines = split_line_content(line)
                        for split_line in split_lines:
                            if any(protein.lower() in split_line.lower() for protein in self.proteins if not isinstance(protein, tuple)) or \
                               any(p.lower() in split_line.lower() for protein in self.proteins if isinstance(protein, tuple) for p in protein):
                                found_lines.append("CORE PROTEIN: " + split_line)
                            else:
                                found_lines.append(split_line)
                        found_lines.append("")
                        found_lines.append("")

                # Add core region
                found_lines.append("\nCore Protein Region:")
                found_lines.append("")
                for line in lines[region_start_index:region_end_index + 1]:
                    split_lines = split_line_content(line)
                    for split_line in split_lines:
                        if any(protein.lower() in split_line.lower() for protein in self.proteins if not isinstance(protein, tuple)) or \
                           any(p.lower() in split_line.lower() for protein in self.proteins if isinstance(protein, tuple) for p in protein):
                            found_lines.append("CORE PROTEIN: " + split_line)
                        else:
                            found_lines.append(split_line)
                    found_lines.append("")
                    found_lines.append("")

                if following_flanking:
                    found_lines.append("\nFollowing Flanking Genes:")
                    found_lines.append("")
                    for line in following_flanking:
                        split_lines = split_line_content(line)
                        for split_line in split_lines:
                            if any(protein.lower() in split_line.lower() for protein in self.proteins if not isinstance(protein, tuple)) or \
                               any(p.lower() in split_line.lower() for protein in self.proteins if isinstance(protein, tuple) for p in protein):
                                found_lines.append("CORE PROTEIN: " + split_line)
                            else:
                                found_lines.append(split_line)
                        found_lines.append("")
                        found_lines.append("")

            return found_lines

        def find_it_gb_CDS(self, filepath):
            found_lines = []
            all_regions = []
            lines = []

            with open(filepath, 'r') as file:
                records = SeqIO.parse(file, "genbank")
                for record in records:
                    for feature in record.features:
                        if feature.type == "CDS":
                            protein_name = feature.qualifiers.get("product", [""])[0]
                            start = int(feature.location.start)
                            stop = int(feature.location.end)
                            lines.append(f"{protein_name}|{start}|{stop}")

            def extract_proteins(line):
                line_lower = line.lower()
                proteins_found = set()
                for protein in self.proteins:
                    if isinstance(protein, tuple):
                        if any(p.lower() in line_lower for p in protein):
                            proteins_found.add(protein[0].lower())
                    else:
                        if protein.lower() in line_lower:
                            proteins_found.add(protein.lower())
                return proteins_found

            def contains_forbidden(line):
                line_lower = line.lower()
                for protein in self.forbidden:
                    if protein.lower() in line_lower:
                        return True
                return False

            # Build index mapping line content → index
            line_to_index = {line.strip(): idx for idx, line in enumerate(lines)}

            # Collect all potential regions
            for i, line in enumerate(lines):
                line_parts = line.strip().split('|')
                if len(line_parts) < 3:
                    continue
                start_position = int(line_parts[1])

                group_line = []
                protein_counter = set()
                min_stop_position = None

                for linex in lines[i:]:
                    linex_parts = linex.strip().split('|')
                    if len(linex_parts) < 3:
                        continue
                    stop_position = int(linex_parts[2])

                    if min_stop_position is None:
                        min_stop_position = stop_position

                    if (stop_position - start_position <= self.length):
                        if contains_forbidden(linex):
                            break
                        proteins_found = extract_proteins(linex)
                        if proteins_found:
                            group_line.append(linex.strip())
                            protein_counter.update(proteins_found)
                            min_stop_position = stop_position
                    else:
                        break

                if len(protein_counter) >= self.number and len(group_line) > 0:
                    all_regions.append((set(group_line), i, len(group_line), start_position, min_stop_position))

            # Sort regions by number of lines descending
            all_regions.sort(key=lambda x: x[2], reverse=True)

            # Ensure regions are unique
            unique_regions = []
            seen_regions = set()

            def is_subset(region, existing_regions):
                region_lines = set(region[0])
                for existing in existing_regions:
                    if region_lines.issubset(existing[0]):
                        return True
                return False

            for region in all_regions:
                if not is_subset(region, unique_regions):
                    region_lines = set(region[0])
                    seen_regions.add(frozenset(region_lines))
                    unique_regions.append(region)

            # Add auxiliary functions
            def get_flanking_lines(region_start_index, region_end_index):
                previous_flanking = []
                following_flanking = []

                for idx in range(max(0, region_start_index - 5), region_start_index):
                    previous_flanking.append(lines[idx].strip())
                
                for idx in range(region_end_index + 1, min(len(lines), region_end_index + 6)):
                    following_flanking.append(lines[idx].strip())
                
                return previous_flanking, following_flanking

            def split_line_content(line):
                parts = line.split('\t', 1)
                if len(parts) > 1:
                    return [parts[0].strip(), parts[1].strip()]
                return [line.strip()]

            def check_specific_distances(region_lines):
                positions = {protein: None for protein_pair in self.specific_distances for protein in protein_pair}
                for (protein1, protein2), max_distance in self.specific_distances.items():
                    found_protein1 = found_protein2 = False
                    for line in region_lines:
                        if isinstance(protein1, tuple):
                            if any(p.lower() in line.lower() for p in protein1):
                                parts = line.split('|')
                                positions[protein1] = (int(parts[1]), int(parts[2]))
                                found_protein1 = True
                        else:
                            if protein1.lower() in line.lower():
                                parts = line.split('|')
                                positions[protein1] = (int(parts[1]), int(parts[2]))
                                found_protein1 = True

                        if isinstance(protein2, tuple):
                            if any(p.lower() in line.lower() for p in protein2):
                                parts = line.split('|')
                                positions[protein2] = (int(parts[1]), int(parts[2]))
                                found_protein2 = True
                        else:
                            if protein2.lower() in line.lower():
                                parts = line.split('|')
                                positions[protein2] = (int(parts[1]), int(parts[2]))
                                found_protein2 = True

                        if found_protein1 and found_protein2:
                            pos1_end = positions[protein1][1]
                            pos2_start = positions[protein2][0]
                            distance = abs(pos2_start - pos1_end)
                            if distance > max_distance:
                                return False
                            positions[protein1] = positions[protein2] = None
                            found_protein1 = found_protein2 = False
                            continue

                return True

            found_lines = []
            for idx, (region_lines, _, _, region_start_position, region_end_position) in enumerate(unique_regions):
                if not check_specific_distances(region_lines):
                    continue

                group_name = f"{self.title} REGION {idx + 1}"
                found_lines.append(f"\n{group_name}:")
                found_lines.append("=" * len(group_name))

                region_start_index = min(line_to_index[line] for line in region_lines if line in line_to_index)
                region_end_index = max(line_to_index[line] for line in region_lines if line in line_to_index)

                previous_flanking, following_flanking = get_flanking_lines(region_start_index, region_end_index)

                if self.check_forbidden_flanks:
                    if any(contains_forbidden(line) for line in previous_flanking + following_flanking):
                        continue

                if previous_flanking:
                    found_lines.append("\nFlanking Genes:")
                    found_lines.append("")
                    for line in previous_flanking:
                        split_lines = split_line_content(line)
                        for split_line in split_lines:
                            if any(protein.lower() in split_line.lower() for protein in self.proteins if not isinstance(protein, tuple)) or \
                               any(p.lower() in split_line.lower() for protein in self.proteins if isinstance(protein, tuple) for p in protein):
                                found_lines.append("CORE PROTEIN: " + split_line)
                            else:
                                found_lines.append(split_line)
                        found_lines.append("")
                        found_lines.append("")

                found_lines.append("\nCore Protein Region:")
                found_lines.append("")
                for line in lines[region_start_index:region_end_index + 1]:
                    split_lines = split_line_content(line)
                    for split_line in split_lines:
                        if any(protein.lower() in split_line.lower() for protein in self.proteins if not isinstance(protein, tuple)) or \
                           any(p.lower() in split_line.lower() for protein in self.proteins if isinstance(protein, tuple) for p in protein):
                            found_lines.append("CORE PROTEIN: " + split_line)
                        else:
                            found_lines.append(split_line)
                    found_lines.append("")
                    found_lines.append("")

                if following_flanking:
                    found_lines.append("\nFollowing Flanking Genes:")
                    found_lines.append("")
                    for line in following_flanking:
                        split_lines = split_line_content(line)
                        for split_line in split_lines:
                            if any(protein.lower() in split_line.lower() for protein in self.proteins if not isinstance(protein, tuple)) or \
                               any(p.lower() in split_line.lower() for protein in self.proteins if isinstance(protein, tuple) for p in protein):
                                found_lines.append("CORE PROTEIN: " + split_line)
                            else:
                                found_lines.append(split_line)
                        found_lines.append("")
                        found_lines.append("")

            return found_lines

        def __str__(self):
            return f"Family: {self.title}, Length: {self.length}, Proteins: {self.proteins}, Number: {self.number}, Forbidden: {self.forbidden}"
