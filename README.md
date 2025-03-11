DNA to Protein Analysis Tool

## Project Overview
This project consists of two main components:
1. **DNA to Protein Parser**: A Python tool for detecting genes in a DNA sequence, translating them into protein sequences, and saving the results into a file.
2. **BLAST Analysis**: The tool uses the BLAST (Basic Local Alignment Search Tool) to compare protein sequences against the NCBI database, find similar sequences, and save the results to a file.

## Features
### DNA to Protein Parser:
- **Gene Detection**: Identifies genes in a DNA sequence by searching for start and stop codons.
- **Protein Translation**: Translates DNA sequences into protein sequences using the genetic code (codon table).
- **File Saving**: Filters and saves protein sequences (with a minimum length of 50 amino acids) to a file.

### BLAST Analysis:
- **Protein Sequence Analysis**: Compares protein sequences against the NCBI database using BLAST.
- **Filtering**: Filters results based on sequence identity (minimum 85% identity).
- **Results Export**: Saves BLAST results, including sequence similarity, alignment length, and E-value, to a text file.

## Technologies Used
- Python 3.x
- Biopython (for BLAST and sequence handling)
- Tkinter (for GUI file selection)
- NCBI BLAST API

## Prerequisites
- Python 3.x installed on your machine.
- Required Python libraries:
  - `requests`
  - `biopython`
  - `tkinter`

You can install the required libraries with the following command:
```
pip install requests biopython
```

## How to Use

### Step 1: Run DNA to Protein Parser
1. **Start the Script**: Run the `protein_parser.py` script.
2. **Select DNA File**: A file dialog will prompt you to select a text file containing the DNA sequence.
3. **Output**: The tool will extract genes, translate them into proteins, and save proteins with a length greater than or equal to 50 amino acids in a file called `proteins_filtered.txt`.
![alt text](https://github.com/AkimzhanRakhimov/ComparingProtein/blob/main/1.png)
### Step 2: Run BLAST Analysis
1. **Start the Script**: Run the `blast_analysis.py` script.
2. **Input**: Enter the number of protein sequences to analyze (from the `proteins_filtered.txt` file generated in Step 1).
3. **Output**: The script will run BLAST against the NCBI database for the selected protein sequences, filter results with at least 85% similarity, and save the results to a file called `blast_results.txt`.

## File Descriptions
- **protein_parser.py**: Contains the logic for gene detection, protein translation, and saving the results.
- **blast_analysis.py**: Contains the logic for running BLAST on protein sequences and saving the results.
- **proteins_filtered.txt**: A text file containing protein sequences filtered by length (>= 50 amino acids).
- **blast_results.txt**: A text file containing the results of the BLAST analysis.

## Example Workflow
1. **Run** `protein_parser.py` to generate the protein sequences from your DNA sequence.
2. **Run** `blast_analysis.py` to compare the proteins against the NCBI database.


