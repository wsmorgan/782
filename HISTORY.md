# Revision History for "782"

## Revision 0.0.7

- Fixed the eigenvectors so that they are properly transposed.

## Revision 0.0.6

- Added several plotting options for the reproduction of the paper.
- Added an outputs/ folder for storing outputs that we want to keep track of.

## Revision 0.0.5

- Fixed another issue in hamiltonian.py to generalize the construction
  of the rectangle integration method to more potentials.
- Added new tests of solve.py and hamiltonian.py.
- Changed the error message in solve.py to a KeyError.
- Removed the redundant default output file name in solve.py's _solve_system.
- Increased number of devisions used for scanning the potential.

## Revision 0.0.4

- Fixed a number of bugs in hamiltonian.py that prevented it from running.
- Updated solve.py to create and solve the system using the Hamiltonian class.

## Revision 0.0.3

- Added the hamiltonian.py file that defines a Hamiltonian class and
  finds its eigenvalues and eigenvectors.
- Removed the __mul__ method from potentials.py.

## Revision 0.0.2

- Added the last unit tests for the initial repository. The code is
  now 100% coverd.

## Revision 0.0.1

- Added unit tests for the three basic potentials (bump down, sho, and kronig penny).

## Revision 0.0.0

The initial commit to the repo.

This includes the basis package that has a single module that defines
the potentials from a cfg file.