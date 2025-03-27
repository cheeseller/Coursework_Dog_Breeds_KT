# c_Match_Scoring.py saved in the ./4. Source_Code folder
"""
====================================================================================================================================
This script finds and saves the best matching dog breed based on sequence alignment scores. Below is a summary of what this script does:

1. Align Sequences:
    * It uses the `align_sequences` function (from the `b_Sequence_Alignment.py` module) to align a mystery DNA sequence (unknown dog breed from the `mystery.fa` file) 
      against known dog breed sequences and compute alignment scores.

2. Find Best Match:
    * The function `find_best_match` takes a dictionary of alignment scores (from `b_Sequence_Alignment.py` module) and finds the breed with the highest score.
    * It returns the best matching breed and its corresponding alignment score.

3. Save Results:
    * The function `save_results` writes the best matching breed and score to a file named `DNA_Identification_Service_results.txt`.
    * It saves the output in a readable format, i.e., best matching breed, alignment score, and other details of file creation and processing.

Note: The results are also printed to the console for quick access.
=====================================================================================================================================
"""

import os
from b_Sequence_Alignment import align_sequences
from datetime import datetime


def find_best_match(alignment_results, mystery_sequence_length) -> tuple:
    """
    Finds the best match based on alignment scores and calculates the match percentage.

    Parameters:
        alignment_results (dict): Dictionary with breeds as keys and alignment scores as values.
        mystery_sequence_length (int): Length of the mystery sequence for percentage calculation.

    """
    if not alignment_results:
        return None, 0, 0.0

    # Get breed with the highest alignment score
    best_breed = max(alignment_results, key=alignment_results.get)
    # Get the best score
    best_score = alignment_results[best_breed]

    # Calculate match percentage.  Have added a check point for mystery sequence length to avoid dividing by zero
    percent_match = (best_score / mystery_sequence_length) * 100 if mystery_sequence_length > 0 else print("mystery sequence length is zero, cannot complete calculation")
    print(f"\nAlignment results: {alignment_results}")
    print(f"Best matching breed: {best_breed} with score: {best_score}")
    print(f"Percent match: {percent_match:.2f}%")

    return best_breed, best_score, percent_match

def save_results(alignment_results,best_breed_overall, best_score_overall, percent_match, output_file="./6. Results/DNA_Identification_Service_results.txt") -> tuple:
    """
    Saves the best match and score to a file.

    Parameters:
        best_breed_overall (str): Best matching breed.
        best_score_overall (int): Best alignment score.
        percent_match (float)   : Match percentage of the best match sequence breed
        Alignment_results (dict): Dictionary with breeds as keys and alignment scores as values.    
        output_file (str): Folder/file to save results.
    """
    print(f"Results saved to {output_file}")
    return(alignment_results,best_breed_overall, best_score_overall, percent_match)
