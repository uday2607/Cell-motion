### Code inspired from Bock, M., Tyagi, A.K., Kreft,
### JU. et al. Generalized Voronoi Tessellation as a
### Model of Two-dimensional Cell Tissue Dynamics. Bull.
### Math. Biol. 72, 1696–1731 (2010).
### https://doi.org/10.1007/s11538-009-9498-3

import numpy as np
import numba as nb

def mwvoro(config):

    ### Store the configuration data to the local variables
    cm = config[0]                 # Number of cells
    nol = config[1]                # Maximum number of overlaps

    ### Explaination of the variable naming
    #
    # c Cell
    # p (cell) Pair
    # t (cell) Tripel (pairwise overlapping)
    # a contact arc candidate
    # f free closure arc
    #
    # i,j,k,l  # cell indices
    # m,mm,... # overlapping cell pair indices
    # n,nn,... # overlap triple (or vertex) indices
    # o,oo,... # pseudo vertex indices
    # p,pp,... # Voronoi arc indices
    # q,qq,... # free, marginal arc indices

    # Realting and linking the various involved objects
    # Cells -> Pairs, Triples
    c2a = np.zeros((cm, nol))      # Indices of the Voronoi arcs around cell
    c2p = np.zeros((cm, nol))      # indices of pairs containing cell
    c2am = np.zeros((1, cm))       # Num of Voronoi arcs around a cell
    c2pm = np.zeros((1, cm))       # Num of pairs containing a cell

    # Pair -> Cells, Triples, Vertices
    pm = 0                         # Num of pairs
    p2s = np.zeros((nol*cm, 7))    # Contact Sphere variables
                                   # 1:Rij 2:Mijx 3: Mijy
                                   # 4:
