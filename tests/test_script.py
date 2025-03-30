import pytest 
import sys
import os 
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import numpy as np
from scipy.stats import norm

# Path set up
code_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "../code"))
sys.path.insert(0, code_dir)

from script import (
    read_fasta,
    find_best_match,
    calculate_p_value,
    breed_name_tag, 
    build_phylogenetic_tree
)

test_dog_breed_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "test_dog_breed.fa"))
test_mystery_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "test_mystery.fa"))

# Fixtures
@pytest.fixture
def added_breeds():
    return read_fasta(test_dog_breed_path)

@pytest.fixture
def added_mystery():
    return read_fasta(test_mystery_path)[0]

@pytest.fixture 
def data_alignment(added_breeds, added_mystery):
    return find_best_match(added_breeds, added_mystery)

# Tests
def test_files(added_breeds, added_mystery):
    assert len(added_breeds) == 3   
    assert len(added_mystery.seq) == 60   

def test_breed_descriptions(added_breeds):
    for record in added_breeds:
        assert "[breed=" in record.description

def test_find_best_match_structure(data_alignment, added_breeds):
    match, seq, score, all_scores = data_alignment
    assert isinstance(match, str)
    assert isinstance(score, float)
    assert len(seq) > 0 
    assert "breed=" in match 

def test_p_value(data_alignment):
    score, all_scores = data_alignment[2], data_alignment[3]
    p_val = calculate_p_value(score, all_scores)
    
    if len(set(all_scores)) == 1:
        assert p_val == 1.0
    else:
        assert 0 <= p_val <= 1   

def test_empty_mystery_file(added_breeds):
    empty_mystery = SeqRecord(Seq(""), id="mystery")
    with pytest.raises(ValueError, match="sequence has zero length"):
        find_best_match(added_breeds, empty_mystery)

def test_breed_name_tag(added_breeds):
    description = added_breeds[0].description 
    extracted_breed = breed_name_tag(description)
    expected_breed = "Azores Cattle dog"
    assert extracted_breed == expected_breed, f"Expected '{expected_breed}', but got '{extracted_breed}'"

def test_tree(added_breeds, added_mystery):
    sequences = added_breeds + [added_mystery]
    tree = build_phylogenetic_tree(sequences)
    assert len(list(tree.get_terminals())) == 4   

if __name__ == "__main__":
    pytest.main(["-v"])
