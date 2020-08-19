Phylorgs
========

Copyright © Guillaume Louvel 2016-2020


Suite of tools for genomics & phylogeny, implemented during my PhD.


Related publication:

* Factors influencing the accuracy in dating single gene trees,
  Guillaume Louvel, Hugues Roest Crollius.
  bioRxiv 2020.08.24.264671 (preprint); doi: <https://doi.org/10.1101/2020.08.24.264671>


## Main modules

- `dendro`: manipulating trees (uses Ete3 and LibsDyogen);
- `seqtools`: manipulating sequences/alignments (uses Biopython);
- `datasci`: data analysis and plots;
- `ensembltools`, `genomicustools`: utils specific to the given data source;
- `genchron`: gene dating project, including snakemake pipelines;
- `duprates`: duplication rates project.

## Visualisation tools

- `genetree_drawer.py`: draw _reconciled_ gene trees;
- `seqtools/printal.py`: print alignments to terminal with colors (by codon/nucl/aa);
- `seqtools/plot_al_conservation.py`: plot alignment with column stats, and side tree;
- `datasci.graphs.plottree`: function for plotting trees with matplotlib.
- `dendro/interactree.py`: simple wrapper around `Ete3` show commands;

## Special dependencies

- [LibsDyogen_py3](https://github.com/DyogenIBENS/LibsDyogen_py3)
- [ToolsDyogen_py3](https://github.com/DyogenIBENS/ToolsDyogen_py3)
- [evosite3D](https://github.com/romainstuder/evosite3d)
- [fluidcondor](https://github.com/gullumluvl/fluidcondor)

Optional dependency: [atavistic-doc-tools](https://gitlab.com/GullumLuvl/atavistic-doc-tools)

## Required executables

### for `genchron`

- [PAML 4 codeml](http://abacus.gene.ucl.ac.uk/software/paml.html)
- [Newick utilities](http://cegg.unige.ch/newick_utils)
- [Beast2](http://www.beast2.org/)

### for `duprates`

- [phylobayes 4.1](https://github.com/bayesiancook/phylobayes) or [phylobayes MPI](https://github.com/bayesiancook/pbmpi)
- [ALE](https://github.com/ssolo/ALE)
- [Generax](https://github.com/BenoitMorel/GeneRax)
