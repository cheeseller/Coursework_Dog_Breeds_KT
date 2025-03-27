# main_Sequence_Service.py saved in the ./4. Source_Code folder

"""
==============================================================================================================================
this is the script that runs all functions related to the sequences service for delivering the outputs as per course quidelines
==============================================================================================================================
"""
import os
# below requires the installation of psutil.  Run pip install psutil from terminal

# Import the developed functions
from a_Input_from_files import read_fasta
from b_Sequence_Alignment import align_sequences
from c_Match_Scoring import find_best_match, save_results
from datetime import datetime

from d_Phylogeny_tree import build_phylo_tree, draw_phylo_tree
import matplotlib.pyplot as plt
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio import Phylo
 

# File paths
#dog_breeds_file = "Dog_Breeds_Project/3. Data/dog_breeds.fa"
#mystery_file = "Dog_Breeds_Project/3. Data/mystery.fa"

def main():
    """
    Runs all scripts:
    a_Input_from_files.py
    b_Sequence_Alignment.py
    c_Match_Scoring.py
    d_Phylogeny_tree.py
    """

    # Define path/file of dog breed and mystery sequences
    dog_breeds_file = "Dog_Breeds_Project/3. Data/dog_breeds.fa"
    mystery_file = "Dog_Breeds_Project/3. Data/mystery.fa"
 
    '''
    ======================================================================================================
    a_Input_from_files - 
    Call defined function to read the fasta files with the dog_breeds and the mystery sequences
    ======================================================================================================
    '''
    # Call the read_fasta function to read the dog breed and mystery sequences
    dog_breeds_by_breed = read_fasta(dog_breeds_file)
    mystery_sequences = read_fasta(mystery_file)
    
    # Print a summary from reading the sequences files
    total_dog_seqs = sum(len(seq_list) for seq_list in dog_breeds_by_breed.values())
    print(f"Total dog breed sequences: {total_dog_seqs}")
    print(f"Loaded mystery sequences: {len(mystery_sequences)}")
 
    '''
    ======================================================================================================
    b_Sequence_Alignment - 
    Call defined function to 
    - align the dog breed sequences against the mystery sequence, 
    - find breed with the best breed match and its associated scoring, 
    - calculate and provide best match percentage
    ======================================================================================================
    '''
    # Below accesses the mystery sequence list and uses the mystery DNA sequence for the alignment process.
    mystery_key = list(mystery_sequences.keys())[0]
    mystery_sequence = mystery_sequences[mystery_key][0]
    
    # Call the align_sequences function from the b_Sequence_Alignment.py module that processes the alignment and scoring
    alignment_results, best_breed_overall, best_score_overall = align_sequences(mystery_sequence, dog_breeds_by_breed)
    
    """
    ======================================================================================================
    c_Match_Scoring -
    This script finds and saves the best matching dog breed based on sequence alignment scores. Below is a summary of what this script does:

    1. Find Best Match  - Takes a dictionary of alignment scores (from `b_Sequence_Alignment.py` module) and finds the breed with the highest score.
    2. Match percentage - Calculate match percentage based on the best breed match score over the length of the mystery sequence.
    2. Save Results     - Writes the best matching breed, score and matching percentage in a txt file `DNA_Identification_Service_results.txt`.
    ======================================================================================================
    """
    # provide the best breed match and its score, calculate percent match based on the best breed match score over the length of the mystery sequence
    percent_match = (best_score_overall / len(mystery_sequence)) * 100
    print(f"\nAlignment results: {alignment_results}\n")
    print(f"Best matching breed: {best_breed_overall} with score: {best_score_overall}")
    print(f"Percent match: {percent_match:.2f}%")
    

    # Define the full file path
    file_path = os.path.join("Dog_Breeds_Project/6. Results/DNA_Identification_Service_results.txt")
    
    
    # Get current date and time
    current_time = datetime.now().strftime("%d-%m-%Y %H:%M")
    save_results(alignment_results, best_breed_overall, best_score_overall, percent_match)
    with open(file_path, "w") as file:
        # Format information that will be written in the text output file    
        #file = open(file_path, "w")
        file.write("\nCourse Work - Dog Breeds project: DNA Identification Service\n")
        file.write("Prepared by: Konstantinos Tyropolis, MSc Bioinformatics (year 1)\n")
        file.write("Results File: DNA_Identification_Service_results.txt, saved in the Dog_Breeds_Project/6. Results folder\n")
        file.write(f"Date and time file created: {current_time}\n")
        file.write("\n=== Sequences Loaded ===\n")
        file.write(f"Total dog breed sequences loaded: {total_dog_seqs} \n")
        file.write(f"Loaded mystery sequences_by_breed: {len(mystery_sequences)}\n")        
        file.write("\n=== Sequence Match Results ===\n")
        file.write(f"Best Breed Match: {best_breed_overall}\n")
        file.write(f"Match Score: {best_score_overall:.2f}\n")
        file.write(f"Match Percentage: {percent_match:.2f}%\n")
        file.write("\n=== Alignment Results ===\n")
      
        #file = open("DNA_Identification_Service_results.txt", "w")
        file.write("The alignment results 'breed and score' are:\n")

        # Convert dictionary items to a list of formatted "breed: score" strings and print in multiple lines
        entries = [f"{breed}: {score}" for breed, score in alignment_results.items()]
        # Write multiple breeds with their score in multiple lines
        for i in range(0, len(entries), 6):
            file.write(", ".join(entries[i:i+6]) + "\n")  # Write 6 entries per line
                
        #file.write(f"The alignment results 'breed and score' are: \n{alignment_results}")
        file.write("\n=== Phylogeny Tree ===\n")
        file.write(f"Image created and saved in the Dog_Breeds_Project/6. Results folder\n")
        file.write(
            "\n"
            "\n"
            "\n"
        )
        file.write(
            "\nThis output is the result of the functions in the:\n"
            "1. a_Input_from_file.py    ---> to read the dog_breeds.fa and mystery.fa files\n"
            "2. b_Sequence_Alignment.py ---> for aligning sequences using the pairwise2 function and creating a dictionary with sequences and their corresponding scores\n"
            "3. c_Match_Scoring.py      ---> to find the best match, its corresponding score and percentage of match.\n"
            "4. d_Phylogeny_tree.py     ---> to build a phylogenetic tree from the dog breed sequences and includes the mystery sequence.\n"
        )
        file.write("\nFor more details about the project, refer to the README file.\n") 

    
    '''
    ======================================================================================================
    d_phylogeny_tree - 
    Call the defined functions of this script to build a phylogenetic tree from the dog breed sequences and include the mystery sequence
    ======================================================================================================
    '''
    # Build the NJ tree including the mystery sequence and save the drawn image
    tree = build_phylo_tree(dog_breeds_file, mystery_file)

    output_image_path = draw_phylo_tree(tree, output_folder="Dog_Breeds_Project/6. Results")

if __name__ == "__main__":
    main()