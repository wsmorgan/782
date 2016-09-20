#!/usr/bin/python
from basis import msg

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
    "-plot": dict(action="store_true",
                  help=("Plot the solution.")),
    "-potential": dict(help=("Path to the file that has the potential parameters.")),
    "-outfile": dict(default="output.dat",
                     help="Override the default output file nome.")
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

    return args

def run(args):
    print("RUNNING",args)

if __name__ == '__main__': # pragma: no cover
    run(_parser_options())
