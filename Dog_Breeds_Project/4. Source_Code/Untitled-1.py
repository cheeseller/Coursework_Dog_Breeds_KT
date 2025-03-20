"""
This function will do the sequence alignment using pairwise2, which is within the BioPython library.
The function will evaluate each sequence read from dog_breeds.fa vs. the mystery sequence (mystery.fa)
and calculate an alignment score. Each time a higher alignment score is calculated, both the alignment 
score and the corresponding sequence (best match) will be saved in a variable.

Includes:
align_sequences: Uses pairwise2.align.globalxx() to align two sequences globally.
find_best_match: Loops through the database and finds the one with the highest alignment score.
"""

"""
Sequence Alignment Script
-------------------------
Aligns a mystery dog sequence against a known database of dog breeds using Biopython's pairwise2.

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
        dog_breed_sequences (dict): {breed_name: sequence}

    Returns:
        tuple: (best_matching_breed, best_alignment_score)
    """
    print("DEBUG: Entered find_best_match function.")  # Debug
    best_matching_breed = None
    best_alignment_score = float('-inf')

    # Show how many breed sequences we have
    print(f"DEBUG: Total dog_breed_sequences to compare: {len(dog_breed_sequences)}")

    for breed_name, breed_sequence in dog_breed_sequences.items():
        print(f"DEBUG: Aligning with breed: {breed_name}...")  # Debug
        alignments = pairwise2.align.globalxx(mystery_sequence, breed_sequence, one_alignment_only=True)
        
        # If no alignments, move on
        if not alignments:
            print(f"DEBUG: No alignments returned for breed: {breed_name}")
            continue
        
        alignment_score = alignments[0][2]
        print(f"DEBUG: Score for {breed_name}: {alignment_score}")  # Debug
        
        if alignment_score > best_alignment_score:
            best_matching_breed = breed_name
            best_alignment_score = alignment_score
            print(f"DEBUG: New best match is {breed_name} with score {alignment_score}")

    print("DEBUG: Finished find_best_match function.")  # Debug
    return best_matching_breed, best_alignment_score
    print(best_matching_breed, best_alignment_score)

if __name__ == "__main__":
    print("DEBUG: Starting Aligned_Sequences script...")  # Debug

    script_directory = os.path.dirname(os.path.abspath(__file__))
    dog_breeds_fasta = os.path.join(script_directory, "../data/dog_breeds.fa")
    mystery_fasta = os.path.join(script_directory, "../data/mystery.fa")

    print(f"DEBUG: dog_breeds_fasta path: {dog_breeds_fasta}")  # Debug
    print(f"DEBUG: mystery_fasta path: {mystery_fasta}")       # Debug

    dog_breed_sequences = read_fasta(dog_breeds_fasta)
    mystery_dog_sequences = read_fasta(mystery_fasta)

    # Print how many sequences were loaded
    print(f"DEBUG: dog_breed_sequences loaded: {len(dog_breed_sequences)}")  # Debug
    print(f"DEBUG: mystery_dog_sequences loaded: {len(mystery_dog_sequences)}")  # Debug

    if not dog_breed_sequences:
        print("Error: Dog breeds FASTA file is empty or missing.")
        exit()

    if not mystery_dog_sequences:
        print("Error: Mystery FASTA file is empty.")
        exit()

    # Show the keys from mystery_dog_sequences
    print(f"DEBUG: mystery_dog_sequences keys: {list(mystery_dog_sequences.keys())}")

    # Extract the first mystery dog's sequence
    mystery_dog_id, mystery_dog_sequence = list(mystery_dog_sequences.items())[0]
    print(f"DEBUG: Mystery dog ID: {mystery_dog_id}, length: {len(mystery_dog_sequence)}")  # Debug

    best_match, best_score = find_best_match(mystery_dog_sequence, dog_breed_sequences)

    if best_match:
        print(f"Best matching breed: {best_match} with alignment score: {best_score}")
    else:
        print("No suitable alignment found.")

    print("DEBUG: End of Aligned_Sequences script.")  # Debug
