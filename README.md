# biggrExtra

__Metabolic pathway modeling based on gene expression__

The package provides an user-fiendly interface between the gene expression data (expression regulation estimates and their errors) and metabolic modeling with tools provided by the BiGG database (http://bigg.ucsd.edu/) and the R package _BiGGR_ (https://www.bioconductor.org/packages/devel/bioc/html/BiGGR.html). The functionalities include SBML model construction, Monte Carlo calculation of reaction regulation estimates and their errors, significance testing, plotting and generation of hypergraphs. Finally, the gene - SBML model data may be fed into a linear inverse model (LIM) to perform a canonical flux balance analysis (FBA): this functionality is however still in development.

__Installation__

```
devtools::install_github('PiotrTymoszuk/biggrExtra')

```

The package requires also a development package _microViz_ (https://github.com/PiotrTymoszuk/microViz). Make sure to install it first.

__Acknowledgements__

Many thanks to the authors of Bioconductor's _BiGGR_ package, Anand K. Gavai and Hannes Hettling. In particular, code of their tool _sbml2hyperdraw()_ was modified to accommodate hypergraphs.
