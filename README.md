
# genelife

Genetic extension to Conway's Game of Life
by John S. McCaskill and Norman H. Packard

Conway's Game of Life (GoL) is extended to include the influence of genes
proliferating as part of the dynamics. The birth of a cell, depends
on the presence of live neighbours and their genes determine the gene
of the newly born cell. The model optionally further breaks the
(semi-totalistic) symmetry of the GoL, in four steps, allowing the genes
to determine rules depending not just on the number of live neighbours.
Genetic and spatial patterns occuring in the simulation are recorded and
activity statistics, genealogies and a portfolio of other graphical
analysis tools provided. 

Further details are available in the publication by the authors:

John S. McCaskill, Norman H. Packard (2019) "Analysing Emergent Dynamics of Evolving Computation in 2D Cellular Automata"
Lect. Notes in Comp. Sci 11934 3-40. (Open Access)
Online version:  https://doi.org/10.1007/978-3-030-34500-6_1
Print version:   https://link.springer.com/content/pdf/10.1007%2F978-3-030-34500-6_1.pdf

Presented as invited talk at:
TPNC: International Conference on Theory and Practice of Natural Computing
Theory and Practice of Natural Computing
8th International Conference, TPNC 2019, Kingston, ON, Canada, December 9–11, 2019, Proceedings
Editors:  Carlos Martín-Vide, Geoffrey PondMiguel, A. Vega-Rodríguez

## Quick introduction to compiling and running the code.

The main version is designed to run from a Jupyter notebook using a python3 and PySDL2 GUI as well as a C library.

Requirements:
-   Python packages numpy, matplotlib, sdl2 (from `pip install PySDL2`)
-   For SDL2: In addition to PySDL2 (pip install) and the SDL2 library, the library SDL2TTF for true type fonts is employed
(see https://www.libsdl.org/projects/SDL_ttf/).
-   The jupyter notebook interface (see https://jupyter.org/install.html).
-   In addition to Jupyter, the notebook extensions nbextensions are used for convenient table of contents and initialization cells; see https://jupyter-contrib-nbextensions.readthedocs.io/en/latest/install.html.

We suggest that users copy the genelife.ipynb and experiment with the program in their own copy.


The C library is compiled on unix simply by typing make, and this also works on OSX (with some differences in the optimization options compared with the Xcode build). The source code for the C library is contained in the directory genelifec. In addition to the main h file there is a short h file genelife_size.h that specifies the lattice size and related constants.  The lattice size can be changed from within the Jupyter notebook genelife.ipynb by changing N in and executing the second cell to launch compilation.  Apart from the help information provided in the notebook, help is obtained by typing h in the graphics window after starting a run, and both the python and C files have extensive documentation on available functions. There are two custom python files for genelife, the main GUI interface genelife.py and the python-C interface file genelife_c_interface.py.

The code has been tested on Mac OSX versions up to Mojave (10.14.6) compiled using Xcode 11.2.1 (Xcode project is `genelife.xcodeproj`) and on linux (e.g. Ubuntu, CentOS) compiled with make.

## Standalone C-only version using xterm graphics (with restricted functionality)

For those who do not have access to python and SDL2, a much more limited standalone terminal application (using only C) has also
been provided, via the single main C file in the genelifecx directory. The Xcode project file will also make this program for OSX. 
On linux, execute the following line to build
```
	make -f Makefile_xterm
```
This produces an executable `build/genelife_xterm` that may be run at the command line.

On most xterm terminals (256 colors), before running, it makes sense to:
(i)  Reduce the lattice size N to 128 (set log2N to 7 in `genelifec/genelife_size.h` before making)
(ii) In the terminal preferences choose a small font (e.g. 10 pt) 
     and adjust the line and character spacing to 0.5 and 1 for square output.
     
The main C file (`genelifecx/genelife_nopython.c`) runs one example for a certain number of steps, and may be edited to run different options.
However, most of the graphical analysis and program options are only available comfortably from the python-SDL2 GUI above.

