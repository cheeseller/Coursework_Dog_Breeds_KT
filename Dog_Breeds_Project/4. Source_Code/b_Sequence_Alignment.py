# b_Sequence_Alignment.py saved in the ./4. Source_Code folder
 
"""
===================================================================================================================
This function will do the sequence alignment using pairwise2, which is within the BioPython library.
The function will evaluate each sequence read from dog_breeds.fa vs. the mystery sequence (mystery.fa)
and calculate an alignment score. Each time a higher alignment score is calculated, both the alignment
score and the corresponding sequence (best match) will be saved in a variable. best match percentage is also calculated

Includes:
align_sequences: Uses pairwise2.align.globalxx() (alternatively .localxx and localms are used) to align two sequences.
find_best_match: Loops through the database and finds the one with the highest alignment score.
stores best match breed and its corresponding alignment score.
=========================================================================================================================
"""

from Bio import pairwise2
from datetime import datetime
 
def align_sequences(mystery_sequence, dog_breed_groups) -> tuple:
    """
    Aligns the full mystery sequence against the sequences from the first two breed groups
    and returns the alignment score for each group, along with the best overall match.
 
    Parameters:
        mystery_sequence (str): Full DNA sequence of the unknown dog.
        dog_breed_groups (dict): Dictionary { breed_name: [sequence1, sequence2, ...] }.
   
    Returns:
        tuple: (alignment_results, best_breed, best_score)
               alignment_results is a dict with entries for each processed breed.
    """
    # Dictionary to store the alignment results
    alignment_results = {} 
    # Counter to track the number of alignments
    count = 1
    #Variable to store the best alignemnt score, at the end of the process this will be the best score overall
    best_score_overall = 0
    #Variable to store the best breed
    best_breed_overall = None

    print("\n ----------------- START OF DNA IDENTIFICATION SERVICE ------------------")

    '''
    Below is the main script for performing sequence alignment and for finding the best sequence match to the sequence in the mystery fasta file.  In summary:     

    * Extracts breed and breed_sequence.
    * Performs sequence alignment using pairwise2 and records the processing time as information only.  
      Note: All three methods of pairwise2 (globalxx, localxx, localms) were used during the code developent and testing.  The final version of the script is using the localxx method as it appears to be faster than the other two which are commended out.
    * Calculates the alignment score and stores it.
    * Compares each alignment score and provides the best asignment breed and its score.
    ''' 
       
    # for each breed iterate through the dictionary and perform the sequence alignment. Start and end times are used to calculate the processing time of each alignment
    start_time_all = datetime.now()
    for breed, sequences in dog_breed_groups.items():
        breed_best_score = 0
        # Align mystery sequence with each sequence for this breed group 
        for seq in sequences:
            # Provides information on the progress of the alignment process
            print(f"* Alignment: {count} for Breed : {breed} - is progressing. Please wait while processing.....")

            # Defines the start time and further below the end time of the alignment process. Used for information purposes only
            start_time = datetime.now()
            #Script now performs sequence alignment using pairwise2 local method. Addional methods are: globalxx, localms
            alignments = pairwise2.align.localxx(mystery_sequence, seq, one_alignment_only=True, score_only=True) 
            end_time = datetime.now()
            alignment_time = (end_time - start_time).total_seconds()
            print(f"   --------> Alignment with pairwise2 completed. Processing time was {alignment_time:.4f} seconds")
            # For each sequence alignment as per above process, create a dictionary to hold the breed and the alignment score
        '''
        Below lines of code, extract the alignment score for each breed and updates the best score every time there is 
        a higher alignment score.
        '''
        if alignments:
            print("   --------> Scoring extraction is progressing.......")
            # Stores the alignment score for each breed
            score = alignments
            # Saves score every time there is a higher score
            if score > breed_best_score:
                breed_best_score = score
            # Update results in dictionary with best breed match and is score. hold info on the best breed name and the corresponding score
            alignment_results[breed] = breed_best_score
            print(f"   --------> For Alignment {count} for Breed : {breed} the score extracted is: {alignment_results[breed]}")
            count += 1
            if breed_best_score > best_score_overall:
                best_score_overall = breed_best_score
                best_breed_overall = breed
                print(f"   --------> Best matching breed is {best_breed_overall} with a score of {best_score_overall}")
    end_time_all = datetime.now()
    processing_time_all = (end_time_all - start_time_all).total_seconds()
     
    print(f"*\nEnd of Alignment process for all breeds. Sequences alignment processed in {processing_time_all:.3f} seconds ({processing_time_all / 60:.3f} minutes).")
    print("\n------------------- END OF DNA IDENTIFICATION SERVICE -------------------")

    return alignment_results, best_breed_overall, best_score_overall