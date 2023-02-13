from rna_secstruct_design.constraints import (
    MaxRepeatingConstraint,
    MaxRepeatingIncreaseConstraint,
    MaxGCStretchConstraint,
    MaxGCStretchIncreaseConstraint,
)


def test_max_repeating_constraint():
    """Test that the max repeating constraint works"""
    max_repeating_constraint = MaxRepeatingConstraint(2)
    assert max_repeating_constraint.satisifes("AACA")
    assert not max_repeating_constraint.satisifes("AAAAA")
    assert max_repeating_constraint.satisifes("AAGGCCUU")


def test_max_repeating_increase_constraint():
    """Test that the max repeating increase constraint works"""
    max_repeating_increase_constraint = MaxRepeatingIncreaseConstraint(2, "GGG")
    assert max_repeating_increase_constraint.satisifes("AAGGCCUU")
    assert max_repeating_increase_constraint.satisifes("AAGGGCCU")


def test_max_gc_stretch_constraint():
    """Test that the max gc stretch constraint works"""
    max_gc_stretch_constraint = MaxGCStretchConstraint(2)
    assert max_gc_stretch_constraint.satisifes("CAGGAAAACCUG", "((((....))))")
    assert max_gc_stretch_constraint.satisifes("GGGGAAAACCCC", "((((....))))") == False


def test_max_gc_stretch_increase_constraint():
    """Test that the max gc stretch increase constraint works"""
    con = MaxGCStretchIncreaseConstraint(2, "GGGAAACCC", "(((...)))")
    assert con.satisifes("CAGGAAAACCUG", "((((....))))")
    assert con.satisifes("GGGGAAAACCCC", "((((....))))") == False
    assert con.satisifes("GGGAAAAAUCCC", "((((....))))")
