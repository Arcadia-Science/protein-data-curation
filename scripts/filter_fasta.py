import re
from collections import defaultdict
import sys

def extract_gene_and_version(header):
# Pattern 1: gene=H3F3A_2
    gene_match = re.search(r'gene=([^_\s]+)_\d+', header)
    if gene_match:
        gene_id = gene_match.group(1)
        version = re.search(r'gene=' + gene_id + r'_(\d+)', header).group(1)
        return gene_id, version, header

    # Pattern 2: [gene=LOC125944225] [protein_id=XP_049520507.1]
    loc_match = re.search(r'\[gene=(LOC\d+)\]', header)
    protein_id_match = re.search(r'\[protein_id=(XP_\d+\.\d+)', header)
    if loc_match and protein_id_match:
        gene_id = loc_match.group(1)
        version = protein_id_match.group(1)
        return gene_id, version, header

    # Pattern 3: TRINITY
    trinity_match = re.search(r'TRINITY_(DN\d+_c\d+_g\d+)(_i\d+)', header)
    if trinity_match:
        gene_id = trinity_match.group(1)
        version = trinity_match.group(2)
        return gene_id, version, header
        
    # Pattern 4: TRANSDECODER
    transdecoder_match = re.search(r'GENE\.([^\s]+)~~([^\s]+\.p\d+)', header)
    if transdecoder_match:
        gene_id = transdecoder_match.group(1)
        version = transdecoder_match.group(2)
        return gene_id, version, header

    return None, None, header

def read_fasta_file(file_name):
    with open(file_name, "r") as f:
        header = None
        sequence = ""
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if header:
                    yield header, sequence
                header = line
                sequence = ""
            else:
                sequence += line
        if header:
            yield header, sequence

def filter_fasta(input_file, output_file):
    gene_dict = defaultdict(list)
    
    for header, sequence in read_fasta_file(input_file):
        gene_id, version, full_header = extract_gene_and_version(header)
        if gene_id and version:
            gene_dict[gene_id].append((version, len(sequence), full_header, sequence))

    to_remove = set()
    
    for gene_id, entries in gene_dict.items():
        sorted_entries = sorted(entries, key=lambda x: x[1], reverse=True)
        for i in range(1, len(sorted_entries)):
            _, _, header_to_remove, _ = sorted_entries[i]
            to_remove.add(header_to_remove)

    with open(output_file, "w") as out:
        for header, sequence in read_fasta_file(input_file):
            if header not in to_remove:
                out.write(header + "\n" + sequence + "\n")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: filter_fasta.py <input_fasta> <output_fasta>")
        sys.exit(1)
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    filter_fasta(input_file, output_file)
