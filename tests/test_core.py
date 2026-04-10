"""Tests for the core bioinformatics modules."""
import pytest
from ppigfinder.core.orf_finder import ORFFinder

def test_orf_finder_initialization():
    """Test if the ORFFinder initializes correctly."""
    finder = ORFFinder(min_aa_length=30)
    assert finder.min_aa_length == 30
