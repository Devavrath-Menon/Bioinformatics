import os
from Bio import SeqIO
from Bio.Align import substitution_matrices

def blosum_score(residue1, residue2, blosum_matrix):
    residue1, residue2 = residue1.upper(), residue2.upper()
    pair = (residue1, residue2)
    return blosum_matrix.get(pair) or blosum_matrix.get(pair[::-1], 0)  # Efficient retrieval with fallback

def process_fasta(fasta_file, motifs, compulsory_matches, blosum_matrix, output_directory, summary_file_path, threshold=0.75):
    records = list(SeqIO.parse(fasta_file, "fasta"))

    # Ensure output directory exists
    os.makedirs(output_directory, exist_ok=True)

    # Use 'with' to handle file opening/closing automatically
    output_files = {}
    for motif_name in motifs:
        output_files[motif_name] = {
            'full': open(os.path.join(output_directory, f"{motif_name}_full_matches.fasta"), 'w'),
            'partial': open(os.path.join(output_directory, f"{motif_name}_partial_matches.fasta"), 'w'),
            'combined': open(os.path.join(output_directory, f"{motif_name}_combined_matches.fasta"), 'w')
        }

    with open(summary_file_path, 'w') as summary_file:
        for record in records:
            if "HAECO" not in record.id:
                continue  # Process only sequences with 'HAECO' in the header

            sequence = str(record.seq).upper()
            summary_file.write(f">{record.id}\n")
            
            for motif_name, position_list in motifs.items():
                compulsory = compulsory_matches.get(motif_name, [])
                window_size = len(position_list)

                for i in range(len(sequence) - window_size + 1):
                    window = sequence[i:i + window_size]

                    # Verify compulsory matches using all()
                    if not all(window[j] in position_list[j] for j in compulsory):
                        continue  # Skip if compulsory match criteria are not met

                    # Calculate match percentage
                    match_count = sum(window[j] in position_list[j] for j in range(window_size))
                    match_percentage = match_count / window_size

                    if match_percentage >= threshold:
                        # Calculate BLOSUM scores only for mismatches
                        mismatch_scores = [
                            blosum_score(window[j], position_list[j][0], blosum_matrix)
                            for j in range(window_size)
                            if window[j] not in position_list[j]
                        ]
                        num_mismatches = len(mismatch_scores)
                        average_blosum_score = sum(mismatch_scores) / num_mismatches if num_mismatches > 0 else 0

                        # Highlight mismatches in output
                        highlight_window = ''.join(
                            f'({char})' if char not in position_list[j] else char
                            for j, char in enumerate(window)
                        )

                        summary_file.write(f"  {motif_name} match at position {i}: {highlight_window} "
                                           f"with {match_percentage * 100:.1f}% match and average BLOSUM score {average_blosum_score:.2f}\n")
                    
                        # Write to appropriate output FASTA files
                        header = f">{record.id}_pos{i}_{match_percentage * 100:.1f}%match\n{highlight_window}\n"
                        output_files[motif_name]['combined'].write(header)
                        if match_percentage == 1.0:
                            output_files[motif_name]['full'].write(header)
                        else:
                            output_files[motif_name]['partial'].write(header)

    # Close all output files automatically with 'with'
    for files in output_files.values():
        for file in files.values():
            file.close()

# Load the BLOSUM matrix
blosum_matrix = substitution_matrices.load("BLOSUM62")

# Define the directory where you want to save the output FASTA files
output_directory = r"C:\Users\Devavrath\Desktop\Bioinformatics Textbooks\Nov -Research Project\March\One Domain"

# Process the FASTA file
process_fasta(fasta_file, motifs, compulsory_matches, blosum_matrix, output_directory, output_file)
