from Bio import SeqIO
from Bio import Align
import numpy as np
from scipy.stats import norm
from Bio import Phylo
import matplotlib.pyplot as plt 

#Function reads sequences from a FASTA file
def read_fasta(file_path):
    
    """
    Reads sequences from a FASTA file as returns a list of sequences
    """
    sequences = list(SeqIO.parse(file_path, "fasta"))
    return sequences 

# Parameters set to BLAST defaults 
aligner = Align.PairwiseAligner() 
aligner.match_score, aligner.mismatch_score, aligner.open_gap_score, aligner.extend_gap_score = 1.0, -2.0, -5, -2

#Function finds the closest match 
def find_best_match(dog_breed_sequences, mystery_sequence):
    
    """
    Aligns and compares the mystery sequence with the dog breed dataset and returns the closest match and score
    """
    highest_score = -1
    closest_match = ""
    closest_sequence = "" # Stores sequence of the closest match 
    
    for record in dog_breed_sequences:
        alignments = aligner.align(mystery_sequence.seq, record.seq)
        score = alignments[0].score 
        
        #
        if score > highest_score:
            highest_score = score 
            closest_match = record.description
            closest_sequence = record.seq
    
    return closest_match, closest_sequence, highest_score
    
#Reads sequences from FASTA files under data folder
file_path1 = "data/dog_breeds.fa"
file_path2 = "data/mystery.fa" 

dog_breed_sequences = read_fasta(file_path1)
mystery_sequence = read_fasta(file_path2)[0]

# Finds the closest match
closest_match, closest_sequence, similarity_score = find_best_match(dog_breed_sequences, mystery_sequence)

# Prints and saves the output to output.txt file in results folder 

output = f"Closest match: {closest_match}\nSequence: {closest_sequence}\nSimilarity score: {similarity_score}\n"
print(output)
with open("results/Output.txt", "w") as f:
    f.write(output)
print(f"Results saved to: results/Output.txt")

# Function calculates p-value

def calculate_p_value(highest_score, all_scores):
    all_scores = np.array(all_scores)
    mean_score = np.mean(all_scores) 
    std_dev =np.std(all_scores)

    if std_dev == 0: 
        return 1.0 if highest_score <= mean_score else 0.0
    
    z_score = (highest_score - mean_score) / std_dev
    return 1 - norm.cdf(z_score)

all_scores = [aligner.align(mystery_sequence.seq, record.seq)[0].score for record in dog_breed_sequences]

closest_match = max(all_scores)

p_value = calculate_p_value(closest_match, all_scores)
print(f"P-value for closest match: {p_value}")

