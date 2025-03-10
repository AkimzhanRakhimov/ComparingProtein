import ssl
from Bio import SeqIO
from Bio.Blast import NCBIWWW, NCBIXML
from io import StringIO

# Disable SSL certificate verification
ssl._create_default_https_context = ssl._create_unverified_context

# Reading protein sequences from a file (FASTA)
def read_proteins_from_file(file_path):
    return [str(record.seq) for record in SeqIO.parse(file_path, "fasta")]

# Remove stop codons from the sequence
def remove_stop_codons(protein_sequence):
    return protein_sequence.replace("*", "")

# Sending BLAST request
def blast_protein_sequence(protein_sequence):
    result_handle = NCBIWWW.qblast("blastp", "nr", protein_sequence)
    return result_handle.read()

# Parsing BLAST results and saving them to a file
def parse_and_save_results(blast_results, output_file):
    blast_record = NCBIXML.read(StringIO(blast_results))
    with open(output_file, "a") as f:
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                identity_percent = (hsp.identities / hsp.align_length) * 100
                if identity_percent >= 85:  # Similarity filter
                    f.write(f"{alignment.title}\t{identity_percent:.2f}%\t{hsp.align_length}\t{hsp.expect}\n")

# Main function
def main():
    proteins = read_proteins_from_file("proteins_filtered.txt")
    num_proteins = int(input(f"How many proteins to test? (1-{len(proteins)}): "))
    
    with open("blast_results.txt", "w") as f:
        f.write("Name\tSimilarity (%)\tAlignment length\tE-value\n")
    
    for i, protein in enumerate(proteins[:num_proteins]):
        protein = remove_stop_codons(protein)
        print(f"\nRunning BLAST for protein {i+1}/{num_proteins}: {protein[:50]}...")
        blast_results = blast_protein_sequence(protein)
        parse_and_save_results(blast_results, "blast_results.txt")
        print("Results saved to blast_results.txt")

if __name__ == "__main__":
    main()
