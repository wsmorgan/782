#!/usr/bin/python
from basis import msg
from basis.hamiltonian import Hamiltonian
import numpy as np
import matplotlib.pyplot as plt

def _solve_system(potcfg, n_basis, n_solutions, xl = None, xr = None, plot_f = None, outfile=None):
    """Solves the system for the given potential and the desired number of
    basis functions. Output is written to file.

    Args:
        potcfg (str): The path to the `pot.cfg` file.
        n_basis (int): The number of basis functions to use in the solution.
        n_solutions (int): The number of solutions to be returned.
        xl (float, optional): The left most edge of the potential if different 
            from that stored in the `pot.cfg` file.
        xr (float, optional): The right most edge of the potential if different 
            from that stored in the `pot.cfg` file.
        plot_f (bool, optional): True if the system is going to be plotted.
        outfile (str, optional): The path to the desired output file.

    Returns:
       Output is saved to a csv file `1D_potential_sol.csv". If plot_f = 
           True a plot window is returned.
    """

    ham = Hamiltonian(potcfg, n_basis, xl, xr)
    eigen_vals = ham.eigenvals[:n_solutions]
    eigen_vecs = ham.eigenvecs[:n_solutions]

    with open(outfile,"w+") as outf:
        outf.write("Eigenval      Eigenvec\n")

        for i in range(n_solutions):
            temp = [str(eigen_vals[i])+"      "]
            for j in range(len(eigen_vecs)):
                temp.append(str(eigen_vecs[i][j]))

            outf.write(" ".join(temp)+"\n")

    L = abs(ham.domain[1] - ham.domain[0])
    if plot_f == "pot": # pragma: no cover
        xs = np.arange(ham.domain[0],ham.domain[1],0.01)
        Vs = list(map(ham.pot,xs))
        plt.plot(xs,Vs)
        plt.show()

    elif plot_f == "waves": # pragma: no cover
        for i in range(3):
            wave = eigen_vecs[i]
            xs = np.arange(ham.domain[0],ham.domain[1],0.01)
            psi_x = []
            for x in xs:
                sin_x = 0
                for n in range(len(wave)):
                    sin_x += wave[n]*np.sqrt(2./L)*np.sin((n+1)*np.pi*x/L)
                psi_x.append(sin_x)
            plt.plot(xs,psi_x)
            
        plt.show()

    elif plot_f =="en": # pragma: no cover
        ens = []
        ks = []
        squ_well = []
        for n in range(35):
            ens.append(eigen_vals[n]/(np.pi**2))
            ks.append(n/L)
            squ_well.append(n*n/(L**2))
            
        plt.plot(ks[:-1],ens[:-1],'ro')
        plt.plot(ks,squ_well)
        plt.xlim((0.,3.))
        plt.ylim((0.,10.))
        plt.show()
    
def examples():

    """Prints examples of using the script to the console using colored output.
    """
    script = "BASIS: 1D quantum potential solver using basis expansion."
    explain = ("For simple 1D potentials such as the infinite square well, " 
               "kronig-penny, ect. This code produces a numerical solution "
               "using a bisis expansion.")
    contents = [(("Solve the potential in `kp.cfg` using 200 basis functions."), 
                 "solve.py 200 -potential kp.cfg",
                 "This saves the solution to the default 'output.dat'."
                 "file in the current directory.")]
    required = ("REQUIRED: potential config file `pot.cfg`.")
    output = ("RETURNS: plot window if `-plot` is specified; solution "
              "output is written to file.")
    details = ("The plotting uses `matplotlib` with the default configured "
               "backend. If you want a different backend, set the rc config "
               "for `matplotlib` using online documentation. However, many "
               "backends don't play well with the animation (depending on OS "
               "type and version, etc., etc.; so use carefully.")
    outputfmt = ("")

    msg.example(script, explain, contents, required, output, outputfmt, details)

script_options = {
    "N": dict(default=100, type=int,
              help=("Specifies the number of basis function to be used.")),
    "-plot": dict(help=("Plot the potential (pot), the wave functions (wave), "
                        "the energies (en).")),
    "-potential": dict(help=("Path to the file that has the potential parameters.")),
    "-outfile": dict(default="output.dat",
                     help="Override the default output file nome."),
    "-solutions": dict(default = 10, type=int,
                       help="The number of solutions to be written to file."),
    "-left_edge": dict(default = None, type=float,
                       help="Override the left most edge of the potential "
                       "that has diffined in potential file."),
    "-right_edge": dict(default = None, type=float,
                       help="Override the right most edge of the potential "
                       "that has diffined in potential file.")    
    }
"""dict: default command-line arguments and their
    :meth:`argparse.ArgumentParser.add_argument` keyword arguments.
"""

def _parser_options():
    """Parses the options and arguments from the command line."""
    #We have two options: get some of the details from the config file,
    import argparse
    from basis import base
    pdescr = "1D Quantum Potential Solver."
    parser = argparse.ArgumentParser(parents=[base.bparser], description=pdescr)
    for arg, options in script_options.items():
        parser.add_argument(arg, **options)
        
    args = base.exhandler(examples, parser)
    if args is None:
        return

    if args["plot"] != None:
        args["plot"] = args["plot"].lower()
    
    return args

def run(args):

    if not args["potential"]:
        raise KeyError("A potential file must be provided using the -potential flag.")

    elif args["plot"]:
        _solve_system(args["potential"], args["N"], args["solutions"], xl=args["left_edge"]
                      ,xr=args["right_edge"], outfile = args["outfile"], plot_f = args["plot"])
    else:
        _solve_system(args["potential"], args["N"], args["solutions"], xl=args["left_edge"]
                      ,xr=args["right_edge"], outfile = args["outfile"])

if __name__ == '__main__': # pragma: no cover
    run(_parser_options())
