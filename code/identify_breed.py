from Bio import SeqIO
from Bio import Align

#Function for reading sequences from a FASTA file
def read_fasta(file_path):
    
    sequences = list(SeqIO.parse(file_path, "fasta"))
    return sequences 

#Function for finding the closest match 
def find_best_match(dog_breed_sequences, mystery_sequence):

    highest_score = -1
    closest_match = ""

    aligner = Align.PairwiseAligner()
    aligner.match_score = 1.0 

    for record in dog_breed_sequences:
        alignments = aligner.align(mystery_sequence.seq, record.seq)
        score = alignments[0].score 
        if score > highest_score:
            highest_score = score 
            closest_match = record.description 
    return closest_match, highest_score
    
#Reading sequences from FASTA files under data folder
file_path1 = "data/dog_breeds.fa"
file_path2 = "data/mystery.fa" 

dog_breed_sequences = read_fasta(file_path1)
mystery_sequence = read_fasta(file_path2)[0]

#Finding the closest match
closest_match, similarity_score = find_best_match(dog_breed_sequences, mystery_sequence)

#Printing the results
print(f"Closest match: {closest_match}")
print(f"Similarity score: {similarity_score}")

#Saving the output to output.py file in results folder 



