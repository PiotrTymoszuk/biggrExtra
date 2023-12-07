# biggrExtra

__Metabolic pathway modeling based on gene expression__

The package provides an user-fiendly interface between the gene expression data (expression regulation estimates and their errors) and metabolic modeling with tools provided by the [BiGG database](http://bigg.ucsd.edu/)[^1] and the R package [BiGGR](https://www.bioconductor.org/packages/devel/bioc/html/BiGGR.html)[^2]. The functionalities include SBML model construction, Monte Carlo calculation of reaction regulation estimates and their errors, significance testing, plotting and generation of hypergraphs. Finally, the gene - SBML model data may be fed into a linear inverse model (LIM) to perform a canonical flux balance analysis (FBA): this functionality is however still in development.

__Installation__

The package requires also a development package _microViz_ (https://github.com/PiotrTymoszuk/microViz). Make sure to install it first.

```
devtools::install_github('PiotrTymoszuk/microViz') ## dependency microViz

devtools::install_github('PiotrTymoszuk/biggrExtra')

```

## Terms of use

The package is available under a [GPL-3 license](https://github.com/PiotrTymoszuk/biggrExtra/blob/main/LICENSE).


## Contact

The package maintainer is [Piotr Tymoszuk](mailto:piotr.s.tymoszuk@gmail.com).


__Acknowledgements__

Many thanks to the authors of Bioconductor's _BiGGR_ package, Anand K. Gavai and Hannes Hettling. In particular, code of their tool `sbml2hyperdraw()` was modified to accommodate hypergraphs.

[^1]: King ZA, Lu J, Dräger A, Miller P, Federowicz S, Lerman JA, Ebrahim A, Palsson BO, Lewis NE. BiGG Models: A platform for integrating, standardizing and sharing genome-scale models. Nucleic Acids Research (2016) 44:D515–D522. doi: 10.1093/NAR/GKV1049

[^2]: Gavai AK, Supandi F, Hettling H, Murrell P, Leunissen JAM, Van Beek JHGM. Using Bioconductor Package BiGGR for Metabolic Flux Estimation Based on Gene Expression Changes in Brain. PLOS ONE (2015) 10:e0119016. doi: 10.1371/JOURNAL.PONE.0119016
