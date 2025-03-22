# c_Match_Scoring.py saved in the ./Course Work/4. Source code folder

"""
This script finds and saves the best matching dog breed based on sequence alignment scores. HBelow is a summary of what this script does:

    1. Align Sequences:
    * It uses the align_sequences function (from the b_Sequence_Alignment.py module) to align a mystery DNA sequence (unknown dog breed from mystery.fa file) 
    against known dog breed sequences and compute alignment scores.

    2. Find Best Match:
    * The function find_best_match takes a dictionary of alignment scores (from b_Sequence_Alignment.py module) and finds the breed with the highest score.
    * It returns the best matching breed and its corresponding alignment score.

    3. Save Results:
    * The function save_results writes the best matching breed and score to a file named DNA_Identification_Service_results.txt
    * It saves the output in a readable format, ie., best Matching Breed, Alignment Score, including other details of file creation and processing

    Note: The results are also printed to the console for quick access.
"""

import os

def find_best_match(alignment_results):
    """
    Finds the best match based on alignment scores.

    Parameters:
        alignment_results (dict): Dictionary with breed names as keys and alignment scores as values.

    Returns:
        tuple: (best_breed, best_score)
    """
    if not alignment_results:
        return None, 0

    best_breed = max(alignment_results, key=alignment_results.get)
    best_score = alignment_results[best_breed]

    return best_breed, best_score

def save_results(best_breed, best_score, output_file="./Dog_Breeds_Project/6. Results/DNA_Identification_Service_results.txt"):
    """
    Saves the best match and score to a file.

    Parameters:
        best_breed (str): Best matching breed name.
        best_score (int): Alignment score.
        output_file (str): folder/File to save results.
    """
    with open(output_file, "w") as file:
        file.write(f"Best Matching Breed: {best_breed}\n")
        file.write(f"Alignment Score: {best_score}\n")
    
    print(f"Results saved to {output_file}")


#Include below is a separate .py script

if __name__ == "__main__":
    from a_Input_from_files import read_fasta
    from b_Sequence_Alignment import align_sequences

    # Load sequences
    dog_breeds_file = "./Dog_Breeds_Project/3. Data/dog_breeds.fa"
    mystery_file = "./Dog_Breeds_Project/3. Data/mystery.fa"

    dog_breeds = read_fasta(dog_breeds_file)
    mystery = read_fasta(mystery_file)

    if mystery:
        mystery_sequence = list(mystery.values())[0]
        alignment_results = align_sequences(mystery_sequence, dog_breeds)

        best_breed, best_score = find_best_match(alignment_results)
        print(f"Best match: {best_breed} with score {best_score}")

        # Save the results
        save_results(best_breed, best_score)
