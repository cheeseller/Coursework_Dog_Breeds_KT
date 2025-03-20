# b_Sequence_Alignment.py saved in the ./Course Work/4. Source code folder
 
"""
This function will do the sequence alignment using pairwise2, which is within the BioPython library.
The function will evaluate each sequence read from dog_breeds.fa vs. the mystery sequence (mystery.fa)
and calculate an alignment score. Each time a higher alignment score is calculated, both the alignment
score and the corresponding sequence (best match) will be saved in a variable.
 
Includes:
align_sequences: Uses pairwise2.align.globalxx() (alternatively .localxx and localms are used) to align two sequences.
find_best_match: Loops through the database and finds the one with the highest alignment score.
stores best match breed and its corresponding alignment score.
 
Note 1: As alignment takes a huge amount of time to complete, I have added steps to:
1. to otimize the code to run only for a few iterations or
2. trim the mystery and breeds sequences to 500bp run the code for a specific breed.
 
above steps were used during teh testing and debuging phases but have been left as comments for future reference/use
 
Note 2: the pairwise2 functions has been used with different paramers (.globalxx, .localxx, .localms) to find out which one is faster in processing
the alignments.  The specific command lines have been left as comments for future reference.  I have used the one that
I beleive is faster than the other two
"""
 
from Bio import pairwise2
from datetime import datetime
 
def align_sequences(mystery_sequence, dog_breed_sequences):
    """
    Aligns the mystery sequence to all known dog breed sequences and returns alignment results.
   
    Parameters:
        mystery_sequence (str): DNA sequence of the unknown dog.
        dog_breed_sequences (dict): Dictionary of {breed_name: sequence}.
   
    Returns:
        dict: Alignment scores for each breed.
    """
   
    alignment_results = {}
    count = 1 # Debug
    best_alignment_score = 0
    print("DEBUG: Start of Alignment") # Debug
    iterations = 0 # Debug
    max_iterations = 3 # Debug
    # Trim the mystery sequence to first 500 bp
    #mystery_sequence = mystery_sequence[:500]
    for breed_name, breed_sequence in dog_breed_sequences.items():
        if iterations >= max_iterations:
            break
        iterations += 1
        print(f"DEBUG: Alignment count is: {count} for Breed Name: {breed_name}") # Debug
        # Trim breed sequence to first 500 bp
        #breed_sequence = breed_sequence[:500]
        start_time = datetime.now()
        alignments = pairwise2.align.globalxx(mystery_sequence, breed_sequence, one_alignment_only=True) # fastest
        #alignments = pairwise2.align.localxx(mystery_sequence, breed_sequence, one_alignment_only=True) # Second fastest
        #alignments = pairwise2.align.localms(mystery_sequence, breed_sequence, 1, 0, 0, 0, one_alignment_only=True) # Third fastest
        end_time = datetime.now()        
        alignment_time = (end_time - start_time).total_seconds()
        print(f"DEBUG: Alignment took {alignment_time:.4f} seconds")  # Debug
 
        if alignments:
            print("DEBUG: Alignments found, scoring extraction is progressing.......") # Debug
            alignment_results[breed_name] = alignments[0][2]  # Extract score
            print(f"DEBUG: For Alignment {count} for Breed Name: {breed_name} the score extracted is: {alignment_results[breed_name]}") # Debug
            alignment_score = alignments[0][2]
            count += 1 # Debug
        else:
            print(f"DEBUG: No alignment for {count} found") # Debug
         
        if alignment_score > best_alignment_score:
            print("DEBUG: Calculating best alignment score.......") # Debug
            best_matching_breed = breed_name
            best_alignment_score = alignment_score
            print(f"DEBUG: Best matching breed is {best_matching_breed} with score {best_alignment_score}") # Debug
    print("DEBUG: End of Alignment") # Debug
    return alignment_results
 

# Will include below is a separate .py script
 
if __name__ == "__main__":
    from a_Input_from_files import read_fasta
 
    # Load sequences
    dog_breeds_file = "./Dog_Breeds_Project/3. Data/dog_breeds.fa"
    mystery_file = "./Dog_Breeds_Project/3. Data/mystery.fa"
 
    dog_breeds = read_fasta(dog_breeds_file)
    mystery = read_fasta(mystery_file)
 
    if mystery:
        mystery_sequence = list(mystery.values())[0]
        results = align_sequences(mystery_sequence, dog_breeds)
        print(f"Alignment results: {results}")


