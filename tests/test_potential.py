"""Tests the evaluation of the potential for single and array-valued argumentss for SHO, bump and KP potentials.
"""

import pytest
from basis.potential import Potential
import numpy as np
import sys

def test_getattr():
    """Tests the attribute re-routing to Potential.params.
    """

    pot = Potential("potentials/kp.cfg")
    assert pot.w == pot.params["w"]

    with pytest.raises(AttributeError):
        pot.dommy
        

def test_KP():
    """Tests the Kronig Penny type potential.
    """

    print(sys.version_info)
    print(sys.version_info > (3,0))
    # if sys.version_info > (3,0):
    #     pot = Potential("potentials/kp_3.cfg")
    # else:
    pot = Potential("potentials/kp.cfg")
        
    params = [(2, 1, -15., 10, 100),
              (1e5, 1e3, -1234., 5, 1e5),
              (1./3, 1./6, -10, 20, 5),
              (np.pi, np.pi/2, -np.sqrt(2), 10, 23)]
    
    for w, s, V0, num, N in params:
        pot.adjust_potential(s=s, w=w, v0=V0, n=num)
        x = w*num-s/2.
        xa = np.linspace(0, w*num, N)
        
        assert pot(x) == 0.
        assert len(pot(xa)) == N
        assert pot(-(5.+num)*w) == 0.
        assert pot(0) == V0
        assert pot(w*num - w +s/2.) == V0
        with pytest.raises(ValueError):
            pot("a")
    pass

def test_sho():
    """Tests the SHO potential.
    """
    pot = Potential("potentials/sho.cfg")
    params = [(2., 1., -15., 100),
              (1e5, 1e3, -1234., 1e5),
              (1./3, 1./6, -10, 5),
              (np.pi, np.pi/2, -np.sqrt(2), 23)]
    
    for a, w, V0, N in params:
        pot.adjust_potential(a=a, shift=w, v0=V0)
        x = w+(a-w)/2
        xa = np.linspace(-a,a,N)
        
        assert pot(x) == V0*(x-w)**2
        assert pot(3./4*a) == V0*(3./4*a -w)**2
        assert len(pot(xa)) == N
        assert pot(-5.*a) == 0
        assert pot(-a) == 0
        assert pot(a) == V0*(a -w)**2
        with pytest.raises(ValueError):
            pot("a")
    
def test_bump():
    """Tests the bump in the square well potential.
    """

    pot = Potential("potentials/bump.cfg")
    params = [(2, 1, -15., 100),
              (1e5, 1e3, -1234., 1e5),
              (1./3, 1./6, -10, 5),
              (np.pi, np.pi/2, -np.sqrt(2), 23)]
    
    for a, w, V0, N in params:
        pot.adjust_potential(a=a, w=w, v0=V0)
        x = w+(a-w)/2
        xa = np.linspace(-a,a,N)
        
        assert pot(x) == 0.
        assert pot(3./4*w) == V0
        assert len(pot(xa)) == N
        assert pot(-5.*a) == 0.
        assert pot(-w) == V0
        assert pot(a) == 0.
        with pytest.raises(ValueError):
            pot("a")
    
            
