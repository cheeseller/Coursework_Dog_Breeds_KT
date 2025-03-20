# Align_sequences.py saved in the 4. Source code folder

"""This function will do the sequence alignment using pairwise2, which is within the  BioPython library. The function will evaluate each sequence read from the dog_breeds.fa file vs. the mystery sequence (mystery.fa) and calculate an alignment score.  Each time a higher alignment score is calculated, both the alignment score and the corresponding sequence (best match) will be saved in a variable.
It includes the following executable steps:
align_sequences: Uses pairwise2.align.globalxx() to align two sequences globally. This function compares two sequences and returns the best alignment.
find_best_match: Loops through the sequences in the database and finds the one with the highest alignment score compared to the test sequence.
"""



"""
Sequence Alignment Script
-------------------------
This script aligns a mystery dog sequence against a database of known dog breeds 
using Biopython's pairwise2.

 Reads FASTA files from Source_Code/data/
Uses pairwise alignment to find the best matching breed
Ensures error handling for missing or empty files
"""
import os
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Input_from_files import read_fasta

def find_best_match(mystery_sequence, dog_breed_sequences):
    """
    Aligns the mystery dog's sequence to all known breeds and finds the best match.

    Parameters:
        mystery_sequence (str): The DNA sequence of the unknown dog.
        dog_breed_sequences (dict): A dictionary of breed sequences {breed_name: sequence}.

    Returns:
        tuple: (best_matching_breed, best_alignment_score)
    """
    best_matching_breed = None
    best_alignment_score = float('-inf')

    for breed_name, breed_sequence in dog_breed_sequences.items():
        alignments = pairwise2.align.globalxx(mystery_sequence, breed_sequence, one_alignment_only=True)
        if not alignments:
            continue
        alignment_score = alignments[0][2]
        if alignment_score > best_alignment_score:
            best_matching_breed = breed_name
            best_alignment_score = alignment_score
            

    return best_matching_breed, best_alignment_score

if __name__ == "__main__":
    script_directory = os.path.dirname(os.path.abspath(__file__))
    dog_breeds_fasta = os.path.join(script_directory, "../data/dog_breeds.fa")
    mystery_fasta = os.path.join(script_directory, "../data/mystery.fa")

    dog_breed_sequences = read_fasta(dog_breeds_fasta)
    mystery_dog_sequences = read_fasta(mystery_fasta)

    if not dog_breed_sequences:
        print("Error: Dog breeds FASTA file is empty or missing.")
        #exit()

    if not mystery_dog_sequences:
        print("Error: Mystery FASTA file is empty.")
        #exit()

    mystery_dog_id, mystery_dog_sequence = list(mystery_dog_sequences.items())[0]

    best_match, best_score = find_best_match(mystery_dog_sequence, dog_breed_sequences)

    if best_match:
        print(f"Best matching breed: {best_match} with alignment score: {best_score}")
    else:
        print("No suitable alignment found.")

