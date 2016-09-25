"""Methods used to setup the Hamiltonian of the system."""

import nupmy as np
from basis import msg
from basis.potetial import Potential

class Hamiltonian(object):
    """Represents the Hamliltonian for a 1D quantum potential.

    Args:
        potcfg (str): path to the potential configuration file. 
        n_basis (int): The number of basis functions to use in the solution.
        xi (float, optional): The left most edge of the potential. If
          not specified then it is assumed to be the left most edge of the
          potential as defined in `potcfg'.
        xf (float, optional): The right most edge of the potential. If
          not specified then it is assumed to be the left most edge of the
          potential as defined in `potcfg'.
    
    Attributes: 
        pot (:obj:`Potential`): A Potential object that
          represents the potential for the 1D quantum system.
        eigenvals (list): The energy eigenvalues for the sysetm.
        eigenvecs (list): The eigenvectors for the system.
        ham (np.ndarray): An array of the hamiltonian.
        domain (list): The region over which the potential is defined.

    Examples:
        >>> from basis.hamiltonian import Hamiltonian
        >>> h = hamiltonian("sho.cfg", 10)
        >>> energy = h.eigenvals()
        >>> eigenvecs = h.eigenves()

    """

    def __init__(self, potcfg, n_basis, xi = None, xf = None):
        self.pot = Potential(potcfg)

        if xi == None:
            xi = _find_xi()
        if xf = None:
            xf = _find_xf():

        self.domain = [xi,xf]
        self.ham = self._construct_ham
        eigensys = np.eigh(self.ham).tolist()
        self.eigenvecs = eigensys[0]
        self.eigenvals = eigensys[1]

    # def __call__(self, value):
    #     """Returns the desired row entries for the hamiltonian.

    def _find_xi(self):
        """Finds the left most edge of the potential.
        """

        xi = None
        for key in self.pot.regions:
            if xi == None:
                xi = min(list(key))
            else:
                xi = min([xi,min(list(key))])
                
        return xi

    def _find_xf(self):
        """Finds the right most edge of the potential.
        """

        xf = None
        for key in self.pot.regions:
            if xf == None:
                xf = max(list(key))
            else:
                xf = max([xf,max(list(key))])

        return xf

    def _construct_ham(self, n_basis):
        """Constructs the hamiltonian matrix for the system.
        
        Args: 
            n_basis (int): The number of basis functions to be used
              in the expansion.
        """
        
        ham = np.array()
        xr, xv, width_b, width_vac = self._find_xrs()
        for n in range(n_basis):
            temp = []
            for m in range(n_basis):
                if n == m:
                    hnm = self._fnn(n, xr, width_b)
                    en = (np.pi**2)*(n**2)/(sum(self.domain)**2)
                    if width_vac != None:
                        hnm += self._fnn(n,xv, width_vac)
                else:
                    hnm = self._fnm(n,m,xr,width_b)
                    en = 0
                    if width_vac != None:
                        hnm += self._fnm(n,m,xv, width_vac)
                temp.append(en+hnm)
            ham.append([temp])

        self.ham = ham
                
    def _find_xrs(self):
        """Finds the mid points of the potential bariers.

        Returns: 
            tuple of lists and ints: The list of the barriers in the
                middle of the well, the list of the barriers 
                that define the vaccum region, the widths of the barriers
                in the well and the width of the vaccum region barriers.
        """

        xr = []
        xv = []
        
        # We need to find the best number of divisions for the system
        # to make sure we aren't missing any bumps. For the average
        # user something like 0.1 will likely suffice, however if the
        # user does something special in their potential then we may
        # need to use a smaller iterative size.
        temp = []
        for key in self.pot.params:
            if key != '__builtins__':
                temp.append(abs(self.pot.params[key]))
        if min(temp) > 1:
            divs = 0.1
        else:
            divs = min(temp)/2.0
            
        xs = np(self.domain[0],self.domain[1]+divs,divs)

        # Now we need to scan through the potential to find the bumps.
        Vt = [self.pot(xi),xi]
        for x in xs:
            if self.pot(x) > Vt[0]:
                Vt = [self.pot(x),x]
            elif self.pot(x) < Vt[0]:
                if 0 in Vt:
                    VL = Vt
                elif (abs(self.domain[0])+abs(self.domain[1])) in Vt:
                    VR = Vt
                else:
                    xr.append(np.mean(Vt[1:]))
                    width_b = abs(Vt[-1]) + divs
            else:
                Vt.append(x)

        # If the ends are widder than the middle then we have a
        # defined vaccum in those regions.
        width_vac = abs(VL[-1]) + divs
        if width_vac != width_b or width_vac != width_b/2.:
            xv.append(VL[-1]-VL[-1]/1000,VR[1])
        elif width_vac != width_b/2.:
            xr.append(np.mean(VL[1:]))
            xr.append(np.mean(VR[1:]))
            width_vac = None
        else:
            xr.append(VL[1])
            xr.append(VR[-1])
            width_vac = None

        return xs, xv, width_b, width_vac
        

    def _fnn(self,n, xr, b):
        """The value of the integral over the basis functions for each x in
        xr for the diagonals of the hamiltonian.
        
        Args:
            n (int): The column and row number for the matrix element.
            xr (list of float): The midpoints of the potential barriers.
            b (float): The width of the potential barrier.

        Returns:
            float: The value of the sum of all hnm from the paper.
        """

        hnm = 0
        L = abs(self.domain[0]) + abs(self.domain[1])
        for x in xr:
            spb = x + b/2.
            smb = x - b/2.
            fnnp = spb/L - np.sin(2*n*np.pi*spb/L)/(2*np.pi*n)
            fnnm = smb/L - np.sin(2*n*np.pi*smb/L)/(2*np.pi*n)
            hnm += self.pot(x)(fnnp-fnnm)

        return hnm

    def _fnm(self, n, m, xr, b):
        """The value of the integral over the basis functions for each x in
        xr for the off diagonals of the hamiltonian.
        
        Args:
            n (int): The row number for the matrix element.
            m (int): The column number for the matrix element.
            xr (list of float): The midpoints of the potential barriers.
            b (float): The width of the potential barrier.

        Returns:
            float: The value of the sum of all hnm from the paper.
        """
        hnm = 0
        L = abs(self.domain[0]) + abs(self.domain[1])
        for x in xr:
            spb = x + b/2.
            smb = x - b/2.
            fnnp = np.sin((m-n)*np.pi*spb*L)/((m-n)*np.pi) - np.sin((n+m)*np.pi*spb/L)/(np.pi*(n+m))
            fnnm = np.sin((m-n)*np.pi*smb*L)/((m-n)*np.pi) - np.sin((n+m)*np.pi*smb/L)/(np.pi*(n+m))
            hnm += self.pot(x)(fnnp-fnnm)

        return hnm
