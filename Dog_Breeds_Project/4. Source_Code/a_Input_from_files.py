# a_Input_from_files.py saved in the ./Course Work/4. Source code folder
'''
Reading the Input Files
This function reads the input FASTA files (dog_breeds.fa and mystery.fa) using BioPython's SeqIO module.
This will load the sequences into Python dictionaries and start processing them when the function is called.
This .py file is located in 4. Source code folder, filename: Input_from_files.py
 
Before calling the function, ensure you have installed bellow libraries via the terminal using pip install"""
'''
 
import os
from Bio import SeqIO
 
def read_fasta(file_path):
    """
    Checks if the file exists and returns an error or an empty dictionary.
    Then reads sequences from a FASTA file and returns sequences in dictionary format.
    Returns a dictionary with Key = sequence ID and Value = sequence.
    Reads a FASTA file and returns a dictionary with sequence IDs as keys and sequences as values.
    """
    if not os.path.exists(file_path):
        print(f"Error: The file '{file_path}' was not found. Please check the path.")
        return {}
 
    sequences = {}
    try:
        for record in SeqIO.parse(file_path, "fasta"):
            sequences[record.id] = str(record.seq)
       
        if not sequences:
            print(f"Error: The file '{file_path}' is empty or not in FASTA format.")
            return {}
 
        return sequences
   
    except Exception as e:
        print(f"Error while reading the file '{file_path}': {e}")
        return {}
 
 
#Will Include below is a separate .py script
 
if __name__ == "__main__":
    # Define file paths
    dog_breeds_file = "./Dog_Breeds_Project/3. Data/dog_breeds.fa"
    mystery_file = "./Dog_Breeds_Project/3. Data/mystery.fa"
 
    # Read sequences from files
    dog_breeds = read_fasta(dog_breeds_file)
    mystery = read_fasta(mystery_file)
 
    # Print summary for validation
    print(f"Loaded {len(dog_breeds)} dog breed sequences.")
    print(f"Loaded {len(mystery)} mystery sequences.")


