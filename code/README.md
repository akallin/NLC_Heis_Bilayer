NLC_Heis_Bilayer
================

Current state of the graph files:
4 lines for each graph.

Line 1)
    a. Graph Identifier (integer # beginning at 0 for the first graph)
    b. Number of sites in the graph
    c. Lattice constant of the graph
    d. Low field (0 for no low field expansion, 1 for low field expansion) [no longer necessary for Heisenberg]

Line 2)
    The realspace coordinates of the sites in the graph.  Pairs of x and y coordinates with no punctuation separating them.

Line 3)
    The graph bonds (adjacency list).  Pairs of sites that are adjacent.  Again with no punctuation.

Line 4)
    Subgraphs.  Pairs of numbers: 
      the first referring to the graph identifier of the subgraph, 
      the second is the number of occurances of that subgraph in the current graph.
