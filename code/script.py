from Bio import SeqIO
from Bio import Align
import numpy as np
from scipy.stats import norm
from Bio import Phylo
import matplotlib.pyplot as plt 
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, DistanceCalculator

# Function to read sequences from a FASTA file
def read_fasta(file_path: str):
    """
    Reads sequences in FASTA format and returns a list of sequences

    Parameters: file_path: path to the FASTA files

    Returns: List of sequences 
    """
    sequences = list(SeqIO.parse(file_path, "fasta"))
    return sequences 

# Pairwise aligner initialized with default BLAST parameters
aligner = Align.PairwiseAligner() 
aligner.match_score = 1.0 
aligner.mismatch_score = -2.0 
aligner.open_gap_score = -5
aligner.extend_gap_score = -2

# Function to find the closest match for the mystery sequence 
def find_best_match(dog_breed_sequences, mystery_sequence) -> tuple[str, str, float, list[float]]:
    """
    Aligns and compares the mystery sequence with the dog breed dataset and returns the closest match and its alignment score.

    Parameters:
    dog_breed_sequences: List of sequences to compare against 
    mystery_sequence: The unknown dog breed sequence 

    Returns: 
    Tuple with: description of closest matching breed, sequence of closest match, highest alignment score, list of all alignment scores 
    """
    highest_score = -1
    closest_match = ""
    closest_sequence = "" # Stores sequence of the closest match 
    all_scores = [] 
    
    for record in dog_breed_sequences:
        alignments = aligner.align(mystery_sequence.seq, record.seq)
        score = alignments[0].score 
        all_scores.append(score)
        
        # Adjust the closest match if a higher score is found 
        if score > highest_score:
            highest_score = score 
            closest_match = record.description
            closest_sequence = record.seq
    
    return closest_match, closest_sequence, highest_score, all_scores

# Function to calculate p-value
def calculate_p_value(highest_score: float, all_scores: list[float]) -> float:
    """
    Calculates the p-value based on the alignment scores 

    Parameters:
    highest_score: The highest alignment score to evaluate 
    all_scores: List of all scores to compare  

    Returns: 
    The p-value of the closest match 
    """
    all_scores = np.array(all_scores)
    mean_score = np.mean(all_scores) 
    std_dev = np.std(all_scores)

    if std_dev == 0: 
        return 1.0 if highest_score <= mean_score else 0.0
    
    z_score = (highest_score - mean_score) / std_dev
    return 1 - norm.cdf(z_score)

# Function to obtain the breed name from description 
def breed_name_tag(description: str) -> str: 
    """
    Obtains the breed name from a sequence description 

    Parameters: 
    description: The sequence description 

    Returns: The breed name obtained from the description 
    """
    if "[breed=" in description:
        return description.split("[breed=")[1].split("]")[0]
    elif "breed " in description:
        return description.split("breed ")[-1].split()[0]
    return "Mystery Sequence"

# Function to build a phylogenetic tree 
def build_phylogenetic_tree(clustered_sequences) -> Phylo.BaseTree:
    """
    Computes a phylogenetic tree based on sequence similarity

    Parameters: clustered_sequences: List of sequences to build the tree

    Returns: A phylogenetic tree 
    """
    alignment = MultipleSeqAlignment([
        SeqRecord(seq=record.seq, id=breed_name_tag(record.description))
        for record in clustered_sequences
    ])

    calculator = DistanceCalculator("identity")
    distance_matrix = calculator.get_distance(alignment)

    tree_constructor = DistanceTreeConstructor()
    tree = tree_constructor.nj(distance_matrix)

    return tree

def main():
    # File paths for the FASTA files
    file_path1 = "data/dog_breeds.fa"
    file_path2 = "data/mystery.fa" 

    # Read sequences from FASTA files 
    dog_breed_sequences = read_fasta(file_path1)
    mystery_sequence = read_fasta(file_path2)[0]

    # Find the closest match
    closest_match, closest_sequence, highest_score, all_scores = find_best_match(dog_breed_sequences, mystery_sequence)

    # Print and save the closest match to output.txt file in results folder 
    output = f"Closest match: {closest_match}\nSequence: {closest_sequence}\nSimilarity score: {highest_score}\n"
    print(output)
    
    # Save the result to output.txt
    with open("results/Output.txt", "w") as f:
        f.write(output)
    print(f"Closest match saved to: results/Output.txt")

    # Calculate p-value
    p_value = calculate_p_value(highest_score, all_scores)
    print(f"P-value for closest match: {p_value}")
    
    # Save p-value to output.txt
    with open("results/Output.txt", "a") as f:
        f.write(f"P-value: {p_value}\n")
    print(f"P-value saved to: results/Output.txt")

    # Group sequences by breed 
    breed_clusters = {}
    for record in dog_breed_sequences:
        breed = breed_name_tag(record.description)
        if breed not in breed_clusters:
            breed_clusters[breed] = []
        breed_clusters[breed].append(record)

    # One sequence per breed 
    clustered_sequences = []
    for breed, sequences in breed_clusters.items():
        breed_ref_seq = sequences[0]
        clustered_sequences.append(breed_ref_seq)

    # Include mystery sequence in the list 
    clustered_sequences.append(mystery_sequence)

    # Compute phylogenetic tree
    tree = build_phylogenetic_tree(clustered_sequences)

    # Plot phylogenetic tree 
    fig, ax = plt.subplots(figsize=(5, 5))
    plt.title("Phylogenetic Tree of Dog Breeds", fontsize=8)

    # Function to change tree labels
    def edited_label_func(clade: Phylo.BaseTree.Clade) -> str:
        """
        Custom label function for the tree

        Parameters: Clade: a node in the tree

        Returns: clade name if it is a terminal node
        """
        return clade.name if clade.is_terminal() else ""

    # Draw tree with custom labels
    Phylo.draw(tree, do_show=False, axes=ax, label_func=edited_label_func)

    # Adjust font size of labels 
    for text in ax.texts:
        text.set_fontsize(7) 

    # Save and display the tree 
    plt.savefig("results/tree.png", format="png")
    plt.show()
    print(f"Phylogenetic tree saved to: results/tree.png")

# Call main function only if this script is run directly
if __name__ == "__main__":
    main()
