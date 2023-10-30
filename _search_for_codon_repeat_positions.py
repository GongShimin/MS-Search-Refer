#Find the exact location of all three and more codon repeats
import os
def find_repeated_codons(sequence, min_repeat=3):
    repeated_codons = []
    current_codon = sequence[0:3]
    repeat_count = 1
    current_position = 0
    for i in range(3, len(sequence), 3):
        next_codon = sequence[i:i+3]
        if next_codon == current_codon:
            repeat_count += 1
        else:
            if repeat_count >= min_repeat:
                repeated_codons.append((current_codon, current_position - 3 * repeat_count + 1, current_position))
            current_codon = next_codon
            repeat_count = 1
        current_position += 3
    if repeat_count >= min_repeat:
        repeated_codons.append((current_codon, current_position - 3 * repeat_count + 1, current_position))
    return repeated_codons

def read_fasta_file(file_path):
    gene_sequences = {}
    current_gene = ""
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith(">"):
                current_gene = line.strip()[1:]
                gene_sequences[current_gene] = ""
            else:
                gene_sequences[current_gene] += line.strip()
    return gene_sequences
gene_sequences = read_fasta_file(os.path.join("CDS_species_reference.fa"))  #Use of FASTA files

fo = open("codon_repeat_species_position.txt", 'w') #Output TXT file of record locations
for gene, sequence in gene_sequences.items():
    sequence = sequence[100:-100]  #Selection of CDS regions based on FASTA files
    repeated_codons = find_repeated_codons(sequence)
    if repeated_codons:
        for codon, start, end in repeated_codons:
            fo.write(f"{gene}\t{codon}\t{start}\t{end+3}\n")

