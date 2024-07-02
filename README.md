# Cones of Weighted and Partial Metrics

*by [Janoš Vidali](https://jaanos.github.io/)*

[![](https://img.shields.io/badge/View_on_GitHub-black?logo=github)](https://github.com/jaanos/cones)

The format for the `.ext` and `.ine` files is described in the [cddlib manual](https://people.inf.ethz.ch/fukudak/cdd_home/cddlibman2021.pdf).

The format of the `.oex` and `.oin` files is based on it, with the difference of the index line, and the last columns representing the incidence number and size of the orbit. Also, the numbers before the index line represent the number of orbits and the number of points, rather than the number of extreme rays/inequalities (which is specified in the `total` line) and columns. In the case of weighted metrics (`(01-)(d)w(wq)met`) the `i,i` columns represent the values *w*<sub>*i*</sub>. When the dimension of the cone is smaller than the dimension of the space it is embedded in, the basis of the orthogonal subspace is given after the `nullspace` keyword.

The `.skg` and `.rdg` files contain two matrices. The first matrix is the list of orbits of edges, represented by a pair of the two extreme rays/inequalities connected by an edge in the orbit. The second matrix is the collapsed adjacency matrix according to the equitable partition of extreme rays or inequalities into orbits, as represented in the `.oex` and `.oin` files. There is no particular order in which the orbits appear in the matrix. In some cases they have been manually marked. For some graphs, there are also pictures of the graph or its complement available.

The `.ead`, `.iad`, `.ecd` and `.icd` files contain a list of extreme rays or inequalities, each followed by a number whose absolute value is the adjacency or incidence number, and a list of adjacent or incident extreme rays/inequalities if the number was positive, or a list of non-adjacent or non-incident extreme rays/inequalities if the number was negative. The extreme rays or inequalities are numbered according to the order in which they appear in the `.ext` and `.ine`
files.

The data was computed using [cddlib](https://people.inf.ethz.ch/fukudak/cdd_home/) by Komei Fukuda, [lrslib](http://www-cgrl.cs.mcgill.ca/~avis/C/lrs.html) by David Avis, and the [polyhedral](https://github.com/MathieuDutSik/polyhedral_common) package for [GAP](http://www.gap-system.org/) by Mathieu Dutour Sikirić, with some additional [functions](ConesOfMetrics.g) for working with cones of metrics.

This is a joint work with Michel Deza and Elena Deza. The article is available at [arXiv](http://arxiv.org/abs/1101.0517).

{% include table.html %}

