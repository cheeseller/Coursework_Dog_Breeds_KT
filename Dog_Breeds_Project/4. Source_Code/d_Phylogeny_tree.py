# d_phylogeny_tree.py saved in the ./4. Source_Code folder
"""
=======================================================================================================================
Phylogeny Builder Script
--------------------------------------------------------
This script builds a phylogenetic tree from the dog breed sequences  and includes the mystery sequence
 
Steps:
1. Read the dog breed FASTA file and group sequences by breed (Read from a_Input_from_files.py).
2. Select one sequence per breed.
3. Read the mystery FASTA file and select the mystery sequence (Read from a_Input_from_files.py).
4. Combine the breed sequences with the mystery sequence (Alignment as per b_Sequence_Alignment.py).
5. Construct a MultipleSeqAlignment object.
6. Compute a distance matrix and build a phylogenetic tree using the Neighbor-Joining algorithm.
7. Print the tree in Newick format and also display it using Matplotlib.
========================================================================================================================
"""
 
from a_Input_from_files import read_fasta
import os
import matplotlib.pyplot as plt
from Bio import Phylo
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor

def build_phylo_tree(dog_fasta, mystery_fasta) -> Phylo.BaseTree.Tree:
    """
    Builds a Neighbor-Joining tree that includes both dog breed sequences and the mystery sequence.
    """
    # Read dog sequences
    breed_groups = read_fasta(dog_fasta)  
    # Select one sequence per breed
    records = [SeqRecord(Seq(seq_list[0]), id=breed) for breed, seq_list in breed_groups.items()]
    # Read mystery sequence
    mystery_groups = read_fasta(mystery_fasta)
    # Select the first sequence from the mystery group
    mystery_seq = list(mystery_groups.values())[0][0]
    # Add mystery sequence
    records.append(SeqRecord(Seq(mystery_seq), id="mystery"))

    # Construct a MultipleSeqAlignment object
    alignment = MultipleSeqAlignment(records)
    # Compute a distance matrix and build a phylogenetic tree using the Neighbor-Joining algorithm
    calculator = DistanceCalculator('identity')
    dm = calculator.get_distance(alignment)
    # Construct a tree using the Neighbor-Joining algorithm
    constructor = DistanceTreeConstructor()
    return constructor.nj(dm)

def draw_phylo_tree(tree, output_folder="./6. Results/") -> None:
    """
    Draws and saves the phylogenetic tree.
    """
    # Create figure
    fig, ax = plt.subplots(figsize=(16, 8))

    # Define label function
    def my_label_func(clade):
        return clade.name if clade.is_terminal() else ""

    # Draw tree
    print("Drawing the phylogeny tree.  Please wait....")
    Phylo.draw(tree, label_func=my_label_func, axes=ax, do_show=False)
    ax.set_title("Dog Breeds Phylogeny Tree")

    # Save image of phylogenetic tree
    output_image_path = os.path.join(output_folder, "phylogenetic_tree.png")
    plt.savefig(output_image_path)
    print(f"Phylogenetic tree image created and saved to: {output_image_path}")

    #plt.show() 
    
    print("\nNeighbor-Joining Phylogeny Tree (ASCII format):\n")
    Phylo.draw_ascii(tree)