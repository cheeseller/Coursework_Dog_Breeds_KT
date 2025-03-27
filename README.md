
README file
  
# Below is a summary of the key project information
a. project scope and the key features to address the requested outcome. The project scope is detailed in the project scope and guidance.doc saved under the 1. Requirements and useful project info folder. 
b. Project scope includes two stretched goals
c. Modules to be used
d. approach and logic used for providing the required output
e. Testing and validation
f. Project structure


# a. Project Scope and the key features to address the requested outcome

## Overview - DNA Identification Service
The DNA Identification Service is designed to process DNA sequences, identify the closest match from a provided database, and output the closest sequence along with its similarity score. The project extends into stretch goals that involve calculating probabilities and building a phylogenetic tree based on the sequence data.

### Key Features:
- **Open and read FASTA files**: using BioPython's SeqIO module the build function loads the sequences of all breeds into a python dictionary that will be processed when the function is called. This .py file is located in 4. Source_Code folder, filename: a_Input_from_files.py.

- **Align breed sequences against a given sequence**: Given a query sequence, the program will compare it against a database of sequences and will provide scoring for each sequence based on sequqnce similarity.  The parwise2 module is used for this purpose. This .py file is located in 4. Source_Code folder, filename: b_Sequence_Alignment.py

- **Best match identification and Output of results**: The function will use info derived from the previous step for calculating the best sequence match and the bet match breed name along with its corresponding score and match percentage. The output is saved in a txt file and includes the closest matching sequence's identifier and the score indicating how closely it matches the input sequence and the percentage of match. It also includes the list of the breed names and their calculated scoring. This .py file is located in 4. Source_Code folder, filename: c_Match_Scoring.py

- **Call and execute the above functions**: The script calls an runs all functions related to the sequences service for delivering the outputs as per course quidelines. This .py file is located in 4. Source_Code folder, filename: EXECUTE - main_Sequence_Service.py

  
# b. The stretched goals this project is trying to address

1. **Probabilities**: Compute the probabilities of the query sequence matching different sequences in the database, and calculate p-values to assess the significance of the matches.
2. **Phylogenetic Tree**: Construct a reconstructed phylogeny based on the sequence alignment results, showing how the sequences are related and where the mystery sequence is located in the tree.

# c. The following modules will be used and need to be installed  via the terminal using the pip install command, copy and paste the following command on terminal:

* pip install biopython pandas, matplotlib, networkx, numpy, openpyxl
* Specifically for the sequence allignment the SeqIO,  pairwise2, datetime, Seq and os modules have been used
* For the stretched goal 2, the MultipleSeqAlignment, SeqRecord, DistanceCalculator, DistanceTreeConstructor and Phylo modules from Biopython have been used

# d. Approach and logic used for providing the required output

## Approach and Logic summary:
The program will follow these steps to derive the solution:
1. Input Validation: Ensure that the input sequence file (`mystery.fa`) and the sequence database (`dog_breeds.fa`) are in the correct format (FASTA).
2. Sequence Alignment: Use the `pairwise2` module to align the test sequence against each sequence in the database, using a localxx alignment strategy.
3. Match Scoring: After performing the alignment, the program will extract the alignment score, which indicates how closely two sequences match providing the calcualted scoring and best match %
4. Closest Sequence Identification: The program will find the sequence in the database with the highest alignment score, indicating the closest match to the input sequence.
5. Output Results: Display and save the results, including the closest match sequence and its alignment score along with other info for the file creation and processing

## Approach and Logic Details.  Basically, below describe the steps taken for going from A --> B
The following steps outline the process flow from reading the input files to generating the desired output. These steps are embedded into the code.

### Step 1: Reading and Validating Input Files (Filename: a_Input_from_files.py)
Objective: Ensure that the `dog_breeds.fa` (database) and `mystery.fa` (query) files are in the correct format and contain valid data.
Details:
 - Read the sequences from both input files using `Bio.SeqIO`.
 - Validate that the files are in FASTA format.
 - Handle any file errors (e.g., file not found, empty file, invalid format).
 -  Validation:
  - Check that the files exist.
  - Ensure that the data in both files is correctly formatted (e.g., valid FASTA file).

### Step 2: Sequence Alignment (Using `pairwise2`)
Objective: Align the test sequence (`mystery.fa`) against the sequences in the database (`dog_breeds.fa`).
Details:
 - Use the `pairwise2` module to perform an alignment of the mystery sequence with each sequence in the dog_breeds FASTA file.database.
 - Store the alignment score and sequence match.
 - The best match will be determined by the highest alignment score.

### Step 3: Match Scoring
Objective: Calculate the similarity score for each alignment.
Details:
 - For each alignment, extract the score and check if it is the best score calculated so far.
 - Track the sequence and breed with the highest score, which will be the closest match.

### Step 4: Output Results
Objective: Output the closest match sequence and the alignment score.
Details:
 - Calculate best breed and highest matching alignment score
 - Calculate the percentage of the match
 - Save the closest matches (breed, score, percentage) in a text file (`DNA_Identification_Service_results.txt`).


### Step 5: (Stretch Goal) Probabilities and p-values
Objective: Calculate the probability of the test sequence matching each sequence in the database.
Details:
 - Pending

### **Step 6: (Stretch Goal) Phylogeny
Objective: Construct a Neighbor-Joining (NJ) phylogenetic tree to analyze the relationship between dog breed sequences and a mystery sequence.
Details:
 - Read FASTA files and group dog sequences by breed, select representative sequences, and integrate the mystery sequence.
 - Using Biopython, perform multiple sequence alignment, calculate pairwise distances, and build an NJ tree that considers sequence similarities. The output is a phylogenetic tree image that is saved under 4. Source code folder, filename: d_phylogeny.py

# e. Testing and Validation (Specific Tests)**

Here are specific tests for each step of the code:

1. **Input Validation Tests**:
   - Test with valid FASTA files containing multiple sequences.
   - Test with an invalid file (e.g., incorrect format or missing file).
   - Test with an empty file.
   
3. **Match Scoring Tests**:
   - Test the match scoring with sequences that are identical, partially matching, or use specific base points of sequences.
   - Ensure the correct sequence is returned with the highest score.
   
4. **Output Tests**:
   - Verify the output file is created correctly.
   - Check that the correct match and score are written to the output file.


# f. Supporting material

- Reading of the course tutorials, slides from lectures about pandas, numpy, networkx and other material from the course
- Under Key information -> coursework data
- https://www.w3schools.com/ that offers explanation and testing of python modules, functions, libraries, python syntax and more
- The following links with useful inforamtion relevant to the objectives of the course work:
      * https://biopython.org/docs/1.75/api/Bio.pairwise2.html
      * https://biopython.org/docs/1.75/api/Bio.Align.html
      * https://biopython.org/docs/latest/Tutorial/chapter_align.html
      * https://www.ncbi.nlm.nih.gov/BLAST/tutorial/Altschul-1.html
      * https://biopython.org/wiki/Phylo
      * https://biopython.org/docs/1.76/api/Bio.Phylo.TreeConstruction.html
      * https://docs.python.org
      * https://genepy.org/ - for doing a lot of tests/exercises to get familiar about commands syntax, etc.
- use of google and AI search to find the correct functions and to mainly get explanations of them

# g. Project structure
To accomodate all the above and for ensuring there is clear documentation, easily accessible and cross referenced, the below project structure has been set up (also included in the project structure.doc saved in the project Structure folder):

# Dog_breeds_project (Parent folder)

## README.md (this file)
## 1.	Requirements and useful project info
      ** Coursework projects.pdf
      ** Coursework_checkpoint.pdf
      ** Coursework_update.pdf
      ** Project scope and guidance.doc
## 2.	Project Structure
      ** Project structure.doc
## 3.	Data
      ** dog_breeds.fa
      ** mystery.fa
## 4.	Source_Code
      ** a_Input_from_files.py
      ** b_Sequence_Alignment.py
      ** c_Match_Scoring.py
      ** d_Phylogeny_tree.py
      ** EXECUTE - main_Sequence_Service.py
## 5.	Tests
      ** test_sequence_matching.py
      ** test_alignment.py
## 6.	Results
      ** DNA_Identification_Service_results.txt
      ** Phylogenetic_tree.png
-----------------------------

