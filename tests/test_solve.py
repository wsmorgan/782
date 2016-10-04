"""Tests the script accuss to the basis solver.
"""
import pytest
import os
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
    
def test_run(capfd):
    """Test that a default solve works properly.
    """

    from basis.solve import run
    f_name ='test_output.dat'
    argv = ["py.test", "2", "-potential", "potentials/bump_2.cfg","-outfile",'test_output.dat',"-solutions","2"]
    args = get_sargs(argv)
    run(args)
    model = open("tests/model_output/inf_square.out","r")
    temp = model.read()
    temp2 = open(f_name,"r").read()
    model.close()
    assert  temp2 == temp
    os.system("rm test_output.dat")

    from basis.solve import run
    f_name ='test_output.dat'
    argv = ["py.test", "2", "-potential", "potentials/bump.cfg","-outfile",'test_output.dat',"-solutions","2"]
    args = get_sargs(argv)
    run(args)
    model = open("tests/model_output/bump_1.out","r")
    temp = model.read()
    temp2 = open(f_name,"r").read()
    model.close()
    assert  temp2 == temp
    os.system("rm test_output.dat")

    # argv = ["py.test", "2", "-potential", "potentials/bump_2.cfg","-outfile",'test_output.dat',"-solutions","2","-plot","pot"]
    # args = get_sargs(argv)
    # run(args)
    # os.system("rm test_output.dat")

    argv = ["py.test", "2"]
    args = get_sargs(argv)
    with pytest.raises(KeyError):
        run(args)    
