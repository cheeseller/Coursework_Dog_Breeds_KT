# a_Input_from_files.py saved in the ./4. Source_Code folder
'''
===================================================================================================================
Reading the Input Files
This function reads the input FASTA files (dog_breeds.fa and mystery.fa) using BioPython's SeqIO module.
This will load the sequences_by_breed into Python dictionaries and start processing them when the function is called.
This .py file is located in 4. Source code folder, filename: a_Input_from_files.py
 
Before calling the function, ensure you have installed bellow libraries via the terminal using pip install
====================================================================================================================
'''
 
import os 
# Regular expression operations module, which will be used for formating the script to identify the breed name from the header of each sequence
import re
# Biopython module for reading and writing sequence files
from Bio import SeqIO 
 
def read_fasta(file_path) -> dict:

    """
    Reads sequences_by_breed from a FASTA file and returns sequences_by_breed in dictionary format.
    Returns a dictionary with Key = breed name (or ID if breed name is not found in header) and Value = sequence.
    """
    if not os.path.exists(file_path):
        print(f"Error: The file '{file_path}' was not found. Please check the path.")
        return {}

    # This is the logic of extracting the breed from the header line. Extracts a string that starts with the argument of '[breed=' and the text after it until the next '[' is found.
    breed_pattern = re.compile(r"\[breed=(.*?)\]") 
    # create an empty dictionary to store the sequences by breed
    sequences_by_breed = {} 
    # create a variable to store the breed name
    breed = None    
    
    """
    Below code reads the fasta files with the SeqIO module and extracts the breed names from the header of the dog_breeds.fa file, 
    their  sequences_by_breed.  If header does no include breed information (as it is the case in the mystery.fa file), the breed ID at the 
    beginning of the sequence line (starting with >) will be used as breed name.
    """
    try:
        with open(file_path) as file: 
            for record in SeqIO.parse(file_path, "fasta"):
                # use of record. attribute to store the full header information
                header = record.description
                # applying the breed extraction logic to get a string related to the breed name.  this will be in the format of 'breed=xxxxx'
                match = breed_pattern.search(header)
                if match:
                    # below is picking up the second part of the string. i.e., excluding the 'breed=', which is then just the breed name
                    breed = match.group(1)
                else:
                    # As header does not include breed information, use the record.id attribute to store the unique identifier (ID) for the sequence start of line after character '>'
                    breed = record.id
            
                # create sequences dictionary to include breed name and sequence
                if breed not in sequences_by_breed:
                    sequences_by_breed[breed] = []
                # append the sequence with the new breed name found
                sequences_by_breed[breed].append(str(record.seq))
            # Provide error message if the sequences_by_breed is empty that might have been caused by a non FASTA format of the dog breeds file. This test is lso part of the test script: Test_a_Input_from_files.py
            if not sequences_by_breed:
                print(f"Error: The file '{file_path}' is empty or not in FASTA format.")
                return {}
  
            return sequences_by_breed
 
    except Exception as e:
        print(f"Error while reading the file '{file_path}': {e}")
        return {}



