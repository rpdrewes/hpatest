# HPA: Hierarchical Partition Analysis tool for Phylogenetic Estimation

HPA is a new technique for estimating a phylogeny (an evolutionary tree) based on input genetic sequence data from a set of taxa. This github contains code for a sample implementation of HPA written in the popular programming language [Python](http://www.python.org).

For comparison, other excellent tools that try to do the same thing include [RAxML](https://sco.h-its.org/exelixis/web/software/raxml/index.html) and [PAUP\*](http://paup.sc.fsu.edu/).

HPA takes input data in the form a a standard [FASTA](https://en.wikipedia.org/wiki/FASTA_format) or similar format, and output phylogenies are produced by HPA in the [Newick](https://en.wikipedia.org/wiki/Newick_format) format. This Newick data can be plotted graphically or analyzed with other means. HPA can estimate branch lengths and do bootstrapping to estimate confidence of its results. The branch length estimation is currently under active development.

HPA is conceptually simple compared to many alternatives and also fairly efficient. Even though this implementation of HPA is in the Python programming language, a language not noted for speed or memory efficiency, it runs fairly quickly and can often do jobs in the same time and memory as alternative C programs while returning comparable results.

## Authors
HPA was conceived and implemented by Guy Hoelzer, Associate Professor of Biology at the University of Nevada, Reno, and his student Rich Drewes.

## Dependencies

* (required) [Dendropy](https://github.com/jeetsukumaran/DendroPy)
* (easily removable) [numpy](http://http://www.numpy.org/)
* (totally optional, but useful) [Matplotlib Python plotting library](https://www.matplotlib.org)
* (optional performance enhancing alternative Python interpreter) [pypy](https://pypy.org)

Dendropy is used for general purpose utility functions such as loading and saving of input and output files and some manipulations of tree structures. Numpy is lightly used for some of the HPA algorithm itself. HPA's dependency on numpy is not very deep though and can be removed if desired with minimal performance cost. Helpful plots can show different aspects of the algorithm, if matplotlib is installed.

Optionally, pypy can be used instead of the standard CPython interpreter to improve performance of HPA. In typical tests on moderate size data sets HPA runs around four times faster with pypy than with CPython, though at the cost of greater memory use.

## Installation

[Coming soon]

## Citing HPA

[Coming soon]
