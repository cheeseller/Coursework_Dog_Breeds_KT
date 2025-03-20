"""Reading the Input Files
We'll start by reading the input FASTA files (dog_breeds.fa and mystery.fa) using BioPython's SeqIO module. This will allow us to load the sequences into Python and start processing them.
File is located in 4.Source code folder , filename: a_Input_from_files.py 
"""

""" Ensure you have installed bellow libraries via the terminal using pip install"""

import os 
import sys
from Bio import SeqIO
#sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
def read_fasta(file_path):
    
    """
    Reads sequences from a FASTA file and returns sequences in dictionary, checks if files exist and returns them or raises error.
    Returns a Dictionary with Sequences ID as key and sequences as value.
    """
    # Check if file exists
    if not os.path.exists(file_path):
        print(f"Error: The file '{file_path}' was not found.")
        return {}
    
    sequences = {}                                            # Create empty dictionary to store sequences
    # Read files and extract sequences using SeqIO module from BioPython
    try:
        for record in SeqIO.parse(file_path, "fasta"):
            sequence_id = record.id                           # Extract sequence ID 
            sequence = str(record.seq)                        # Convert sequence to a string
            sequences[sequence_id] = sequence                 # Store in dictionary
        
        # Check if sequences were read and have content    
        if not sequences:
            print(f"The file '{file_path}' is empty or improperly formatted.")
        return sequences
    
    except Exception as e:
        print(f"Error while reading the file '{file_path}': {e}")
        return {}
    
# Run the function to read sequences from the files    
if __name__ == "__main__":
# Paths to data files
    
    dog_breeds_file = "Dog_Breeds_Project/data/dog_breeds.fa"
    mystery_file = "Dog_Breeds_Project/data/mystery.fa"

    # Read sequences from files
    dog_breeds = read_fasta(dog_breeds_file)
    mystery = read_fasta(mystery_file)
    

    # Print summary of loaded sequences for validation purposes
    print(f"There are {len(dog_breeds)} dog breed sequences.")
    print(f"There are {len(mystery)} mystery dog sequences.")


