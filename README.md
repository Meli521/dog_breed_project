Project: Dog Breed Identification 


GitHub Repository: https://github.com/Meli521/dog_breed_project.git


Project Overview:

This project identifies a dog's breed based on an unknown DNA sequence, using sequence alignment, statistical analysis, and phylogenetic tree construction. 

The process involves:

1. Aligning the the unknown sequence (mystery.fa) against a database of known dog breed sequences (dog_breeds.fa) to find the closest match and similarity score 

2. Calculating the probability (p-value) of the closest match

3. Constructing a phylogenetic tree to visualise the evolutionary relationships between the dog breeds 


Porject Structure:

Dog_breed_project
Data -> contains dog_breeds.fa (list of sequences) and mystery.fa (unknown sequence)
Code -> contains the Python script 
Results -> contains the output
Tests -> contains the tests for the script
Venv -> is the virtual envirment
Requirements.txt -> contains the dependciies needed to run this project 


Input: 
- mystery.fa: file containing the unknown dog breed sequence 
- dog_breeds.fa: file containing the known dog breed sequences 

Example Output:

- Closest match: 
>gb|MW916043.1| [location=mitochondrion] [completeness=complete] [topology=circular] [organism=Canis lupus familiaris] [isolate=eDogPT41] [sub_species=familiaris] [breed=Portuguese Warren dog, small size, smooth hair] [gcode=2] [sex=male] [country=Portugal] Canis lupus familiaris isolate eDogPT41 mitochondrion, complete genome.

- Sequence of the closest match: GTTAATGTAGCTTAATTAATAAAGCAAGGCACTGAAAATGCCAAGATGAGTCGCACGACT
CCATAAACATAAAGGTTTGGTCCTAGCCTTCCTATTAGTTTTTAGTAGACTTACACATGC
AAGCCTCCACGCCC....

- Similarity score: 
17777.0 

- P-value:
0.0299999737188924

- Phylogenetic tree:
saved as 'tree.png' in the results folder 

Instructions:

1. Clone the repository (https://github.com/Meli521/dog_breed_project.git)

2. Set up the Python environment 

3. Install the required libraries found in requirements.txt

4. Ensure the FASTA files 'dog_breeds.fa' and 'mystery.fa' exist in the data/ directory  

5. Run the script in the code/ directory

6. The closest match, the sequence of the closest match, the similarity score, the p-value will be saved in the 'Output.txt' file under results/ directory. The phylogenetic tree will be saved as 'tree.png' in the results/ directory  

External libraries used:
- NumPy
- BioPython
- pytest
- Matplotlib
- SciPy 