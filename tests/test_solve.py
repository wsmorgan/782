"""Tests the script accuss to the basis solver.
"""
import pytest
# The virtual, pseudorandom port is setup as a session fixture in conftest.py
def get_sargs(args):
    """Returns the list of arguments parsed from sys.argv.
    """
    import sys
    sys.argv = args
    from basis.solve import _parser_options
    return _parser_options()    

def test_examples():
    """Makes sure the script examples work properly.
    """
    argv = ["py.test", "-examples"]
    assert get_sargs(argv) is None
    
def test_run():
    """Test that a default solve works properly.
    """

    argv = ["py.test", "200", "-potential", "potentials/kp.cfg"]
    args = get_sargs(argv)
    from basis.solve import run
    assert run(args) == 0

    
