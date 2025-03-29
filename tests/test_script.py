import pytest 
import sys
import os 
from Bio.SeqRecord import SeqRecord
code_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "../code"))
sys.path.insert(0, code_dir)

from script import (
    read_fasta,
    find_best_match,
    calculate_p_value,
    breed_name_tag,
    build_phylogenetic_tree
)

test_dog_breeds_path = os.path.join(os.path.dirname(__file__), "test_dog_breeds.fa")
test_mystery_path = os.path.join(os.path.dirname(__file__), "test_mystery.fa")


def test_functions():
    assert callable(read_fasta)
    assert callable(find_best_match)
    assert callable(calculate_p_value)
    assert callable(breed_name_tag)
    assert callable(build_phylogenetic_tree)

