# biggrExtra

__Prediction of Metabolism Activity with Transcriptome Data__

## Summary

The idea of modeling of activity of metabolism with gene expression information 
is not new. 
Inference of activity of metabolic pathway with RNA sequencing or microarray data 
has a potentially broad application e.g. in cancer biology and immunology studies. 
Availability of simple "gene - reaction association rules" in Recon knowledge models 
of human metabolism [^1] makes it possible to calculate differences of activity 
or activity scores of metabolic pathways with gene expression or differential gene 
expression data. 
Such simple approach, coined "gene - protein - reaction" or GPR, was used by authors 
of the former Bioconductor R package _BiGGR_ to model changes of activity of metabolism 
in brain tissue with transcript levels quantified with a microarray [^2]. 
The GPR strategy assumes that differences in mRNA levels translate directly 
to differences in amounts of their protein products such as enzymes or regulatory 
proteins, which in turn affects activity of metabolic reactions in the most straightforward 
way. 
Of note, the GPR approach ignores manifold biological processes impacting on enzymatic 
activity such as post-translational protein modifications, allosteric regulation, or 
simply stoichometry, i.e. availability of reaction's substrates and products. 
Yet, evaluation of the gene - reaction association rules with the GPR strategy provides 
us with an imperfect but practical screening and hypothesis-generating tool to 
explore the metabolic layer of transcriptome-coded information. 

In our package `biggrExtra`, we implemented the GPR evaluation scheme in a way 
proposed by the R package _BiGGR_ authors [^2] and appended it with statistical inference, 
diagnostic and visualization tools. 
Our package works with manually pre-processed Recon models warranting common formats of gene identifiers 
and gene - reaction association rules delivered as RData objects with the package. 
Currently, two of them, Recon2 [^3] (https://www.ebi.ac.uk/biomodels/MODEL1109130000) by Thiele et al., 
and Recon 2.2 [^4] (https://www.ebi.ac.uk/biomodels/MODEL1603150001) by Swainston et al., 
are available. 
For R code used in harmonization of the reaction annotation data, please visit our public 
GitHub repository [RECON_processing](https://github.com/PiotrTymoszuk/RECON_processing). 

The package enables essentially two types of analyses: 

* __prediction of fold-differences in metabolic reaction activity with differential gene expression data__. 
This task is handled with function `get_regulation()`, which takes a numeric vector of 
differential gene expression estimates and, optionally, a numeric vector of differential 
gene expression errors. If the vector of differential gene expression errors is provided, 
statistical inference for fold-regulation estimates of reaction activity is made, e.g. 
with simulation tests (Monte Carlo)

* __calculation of activity scores of metabolic reactions with mRNA expression levels__ 
with function `get_activity()`. The function takes a numeric matrix or a data frame with 
gene expression information pre-processed by the user (e.g. Z-scores of $log_2$ gene 
expression counts). The activity scores derived from evaluation of the gene - reaction 
association rules for single RNA sequencing of microarray samples can be subsequently 
compared between analysis groups of interest with a range of statistical tests or models. 

Usage examples are provided in below the __Basic usage__ section of the the package's homepage.

## Installation

The easiest way to install the `biggrExtra` package and its 
dependency [fastTest](https://github.com/PiotrTymoszuk/fastTest) used for enrichment 
tests is to use `devtools`: 

```
devtools::install_github('PiotrTymoszuk/fastTest')

devtools::install_github('PiotrTymoszuk/biggrExtra')

```

## Terms of use

The package is available under a [GPL-3 license](https://github.com/PiotrTymoszuk/biggrExtra/blob/main/LICENSE).


## Contact

The package maintainer is [Piotr Tymoszuk](mailto:piotr.s.tymoszuk@gmail.com).


## Acknowledgements

Many thanks to the authors of the former Bioconductor's _BiGGR_ package, Anand K. Gavai and Hannes Hettling. 


## Basic usage 

### Prediction of differences in activity of metabolic reactions with differential gene expression data

#### Analysis data

In this "real-life" analysis example, we will model fold-changes of activity of Recon 2.2. metabolic reactions 
between two molecular subsets of prostate carcinoma differing in expression of collagens and collagen-processing 
genes [^5] with their differential gene expression data. 

In the example, except of our `biggrExtra` package, we'll use the `tidyverse` package bundle for transformation 
of tabular data, plotting, and construction of pipelines. 

For the modeling of reaction's activity with `get_regulation()` function, we'll need three ingredients: 
(1) a vector of estimates of $log_2$ fold-differences in gene expression named with gene's Entrez identifiers (Entrez ID), 
(2) a vector of errors of $log_2$ fold-differences in gene expression named with gene's Entrez identifiers (Entrez ID), and 
(3) a reaction-annotation database of class `reactDB` with ready-to-use gene - reaction association rules as R language objects.

Estimates and errors of $log_2$ fold-differences in gene expression levels between these collagen subsets are 
delivered with our package as `tcga_data` data set. 
The data frame with Recon 2.2 reaction annotation information is provided with our package as `Recon2_2D` object and 
can be easily converted to a `reactDB` database with `as_reactDB()` function.

```r

  # R packages
  
  library(tidyverse)
  library(biggrExtra) ### modeling of metabolism
  
  # analysis data
  
  ## the `tcga_data` data set, delivered with `biggrExtra` package
  ## differential regulation: estimates of log2 fold-changes
  ## with their errors, named with Entrez IDs of the genes

  dge_estimates <- set_names(tcga_data$estimate,
                             tcga_data$entrez_id)

  dge_errors <- set_names(tcga_data$se,
                          tcga_data$entrez_id)

  ## reaction annotation database created from
  ## the `Recon2_2D` data set from `biggrExtra` package

  recon2_2_db <- as_reactDB(Recon2_2D)
  
```
```r
> dge_estimates[1:10]

  106480547   100873672   100873326   106481046   106479697   106480905   106480310   106479086   106479870   106481303 
 0.34097622 -0.07640674  0.22816446  0.04801038  0.11530435  0.38314156  2.13545242  0.00000000  1.72093787 -0.69571057 
 
> dge_errors[1:10]

 106480547  100873672  100873326  106481046  106479697  106480905  106480310  106479086  106479870  106481303 
0.30921224 0.06602058 0.16943733 0.13573124 0.09815172 0.28805994 0.77680608 0.00000000 0.77967968 0.46749132 

```
```r

> recon2_2_db

`reactDB` object with 7785 reactions and 4742 gene - reaction association rules

# A tibble: 7,785 × 6
   id                name                                                    subsystem         gene_association entrez_id exprs     
 * <chr>             <chr>                                                   <chr>             <chr>            <list>    <named li>
 1 R_13DAMPPOX       1,3-Diaminopropane:oxygen oxidoreductase (deaminating)  beta-Alanine met… 314 or 8639 or … <chr [3]> <language>
 2 R_1a_24_25VITD2Hm 1-alpha-Vitamin D-24,25-hydroxylase (D2)                Vitamin D metabo… NA               <NULL>    <NULL>    
 3 R_1a_24_25VITD3Hm 1-alpha-Vitamin D-24,25-hydroxylase (D3)                Vitamin D metabo… NA               <NULL>    <NULL>    
 4 R_1a_25VITD2Hm    1-alpha,24R,25-Vitamin D-hydroxylase (D2)               Vitamin D metabo… NA               <NULL>    <NULL>    
 5 R_1a_25VITD3Hm    1-alpha,24R,25-Vitamin D-hydroxylase (D3)               Vitamin D metabo… NA               <NULL>    <NULL>    
 6 R_1PPDCRp         delta1-piperideine-2-carboxylate reductase, perixosomal Lysine metabolism NA               <NULL>    <NULL>    
 7 R_24_25VITD2Hm    24R-Vitamin D-25-hydroxylase (D2)                       Vitamin D metabo… 1591             <chr [1]> <sym>     
 8 R_24_25VITD3Hm    24R-Vitamin D-25-hydroxylase (D3)                       Vitamin D metabo… 1591             <chr [1]> <sym>     
 9 R_25VITD2Hm       1-alpha-Vitamin D-25-hydroxylase (D2)                   Vitamin D metabo… 1594             <chr [1]> <sym>     
10 R_25VITD3Hm       1-alpha-Vitamin D-25-hydroxylase (D3)                   Vitamin D metabo… 1594             <chr [1]> <sym>     
# ℹ 7,775 more rows
# ℹ Use `print(n = ...)` to see more rows

```

#### Modeling of activity of metabolic reactions

Fold-changes of activity of metabolic reactions in the collagen subsets of prostate cancers 
will be computed with function `get_regulation()` 
with statistical inference (standard errors and bias-corrected and accelerated [BCA] 95% confidence 
intervals) derived from Monte Carlo simulations. 
In the simulations, normal distributions of $log_2$ fold-differences in mRNA levels (distribution mean: 
$log_2$ fold-change estimate, distribution standard deviation: $log_2$ fold-change error) are sampled at random, 
and for each sample the gene - reaction association rules are evaluated with the "gene - protein - reaction" 
GPR strategy.
These task are handled by the code below: 

```r

  ## predictions of fold-changes in reaction activity,
  ## errors derived from Monte Carlo simulations

  reaction_regulation <-
    get_regulation(x = dge_estimates, ### vector of DGE estimates
                   err = dge_errors, ### vector of DGE errors
                   database = recon2_2_db, ### reaction annotation database
                   scale = "log2", ### scale of the DGE estimates
                   err_method = "mc", ### error and p value calculation method
                   return_mc = TRUE, ### reaction activity for each iteration saved
                   n_iter = 2000, ### number of iterations
                   .parallel = TRUE)

```

This step may take several minutes, even if run in parallel (`.parallel = TRUE`). 
The result is an `actiData` object consisting of two elements: 

* `reg`: a data frame with fold-regulation estimates and their inference statistics 
on identity scale, raw p-values and p-values corrected for multiple testing with 
the FDR method (false discovery rate). This data is annotated with reaction names 
and Recon metabolic subsystems

* `mc`: an optional numeric matrix returned only if `return_mc = TRUE` with 
estimates of fold-regulation of reaction activity in particular iterations of the 
simulation algorithm

```r
> reaction_regulation

`actiData` object with activity estimates of 4742 in 94 subsystems

# A tibble: 4,742 × 9
   id             name                                                subsystem fold_reg  error lower_ci upper_ci p_value p_adjusted
   <chr>          <chr>                                               <chr>        <dbl>  <dbl>    <dbl>    <dbl>   <dbl>      <dbl>
 1 R_13DAMPPOX    1,3-Diaminopropane:oxygen oxidoreductase (deaminat… beta-Ala…    1.46  0.130     1.26     1.80   0.0005   0.000797
 2 R_24_25VITD2Hm 24R-Vitamin D-25-hydroxylase (D2)                   Vitamin …    2.37  0.763     1.35     4.54   0.003    0.00438 
 3 R_24_25VITD3Hm 24R-Vitamin D-25-hydroxylase (D3)                   Vitamin …    2.37  0.763     1.35     4.54   0.003    0.00438 
 4 R_25VITD2Hm    1-alpha-Vitamin D-25-hydroxylase (D2)               Vitamin …    0.724 0.0583    0.617    0.850  0.0005   0.000797
 5 R_25VITD3Hm    1-alpha-Vitamin D-25-hydroxylase (D3)               Vitamin …    0.724 0.0583    0.617    0.850  0.0005   0.000797
 6 R_2AMACHYD     2-Aminoacrylate hydrolysis                          Glycine,…    1.91  0.204     1.57     2.37   0.0005   0.000797
 7 R_2AMACSULT    2-Aminoacrylate sulfotransferase                    Methioni…    0.961 0.332     0.536    1.93   0.384    0.412   
 8 R_2HBO         2-Hydroxybutyrate:NAD+ oxidoreductase               Propanoa…    0.975 0.140     0.793    1.36   0.363    0.391   
 9 R_2OXOADOXm    2-Oxoadipate:lipoamde 2-oxidoreductase(decarboxyla… Lysine m…    0.916 0.0225    0.871    0.958  0.0005   0.000797
10 R_34DHOXPEGOX  3,4-Dihydroxyphenylethyleneglycol:NAD+ oxidoreduct… Tyrosine…    2.19  0.520     1.52     3.90   0.0005   0.000797
# ℹ 4,732 more rows
# ℹ Use `print(n = ...)` to see more rows

```

```r
> reaction_regulation$reg[1:10, ]

# A tibble: 10 × 9
   id             name                                                subsystem fold_reg  error lower_ci upper_ci p_value p_adjusted
   <chr>          <chr>                                               <chr>        <dbl>  <dbl>    <dbl>    <dbl>   <dbl>      <dbl>
 1 R_13DAMPPOX    1,3-Diaminopropane:oxygen oxidoreductase (deaminat… beta-Ala…    1.46  0.130     1.26     1.80   0.0005   0.000797
 2 R_24_25VITD2Hm 24R-Vitamin D-25-hydroxylase (D2)                   Vitamin …    2.37  0.763     1.35     4.54   0.003    0.00438 
 3 R_24_25VITD3Hm 24R-Vitamin D-25-hydroxylase (D3)                   Vitamin …    2.37  0.763     1.35     4.54   0.003    0.00438 
 4 R_25VITD2Hm    1-alpha-Vitamin D-25-hydroxylase (D2)               Vitamin …    0.724 0.0583    0.617    0.850  0.0005   0.000797
 5 R_25VITD3Hm    1-alpha-Vitamin D-25-hydroxylase (D3)               Vitamin …    0.724 0.0583    0.617    0.850  0.0005   0.000797
 6 R_2AMACHYD     2-Aminoacrylate hydrolysis                          Glycine,…    1.91  0.204     1.57     2.37   0.0005   0.000797
 7 R_2AMACSULT    2-Aminoacrylate sulfotransferase                    Methioni…    0.961 0.332     0.536    1.93   0.384    0.412   
 8 R_2HBO         2-Hydroxybutyrate:NAD+ oxidoreductase               Propanoa…    0.975 0.140     0.793    1.36   0.363    0.391   
 9 R_2OXOADOXm    2-Oxoadipate:lipoamde 2-oxidoreductase(decarboxyla… Lysine m…    0.916 0.0225    0.871    0.958  0.0005   0.000797
10 R_34DHOXPEGOX  3,4-Dihydroxyphenylethyleneglycol:NAD+ oxidoreduct… Tyrosine…    2.19  0.520     1.52     3.90   0.0005   0.000797

```

```r

> reaction_regulation$mc[1:10, 1:5]

      R_13DAMPPOX R_24_25VITD2Hm R_24_25VITD3Hm R_25VITD2Hm R_25VITD3Hm
 [1,]    1.550390       2.219177       2.219177   0.7603234   0.7603234
 [2,]    1.456067       2.410301       2.410301   0.6971161   0.6971161
 [3,]    1.274681       1.288285       1.288285   0.6991904   0.6991904
 [4,]    1.539576       1.308465       1.308465   0.8823137   0.8823137
 [5,]    1.438337       1.980600       1.980600   0.6767248   0.6767248
 [6,]    1.535606       1.343117       1.343117   0.7506480   0.7506480
 [7,]    1.625061       2.524067       2.524067   0.7819929   0.7819929
 [8,]    1.299729       3.899864       3.899864   0.7096446   0.7096446
 [9,]    1.408779       2.112449       2.112449   0.8277458   0.8277458
[10,]    1.400568       2.953604       2.953604   0.7016552   0.7016552

```

#### $log_2$ fold-regulation estimates of reaction activity and identification of significant effects

Fold-changes of reaction activity on an identity scales are not particularly practical in handling, 
they can be easily transformed with $log_2$ with function `transform_estimates()`. 
This function appends the `reg` data frame of the `actiData` object with new columns 
storing the transformed fold-regulation estimates, errors, and 95% confidence intervals. 
Names of these new columns are preceded with a character string specified by `prefix` argument, 
in our case `"log2_"`. 
The transformation function is provided as `fun` argument: 

```r

  reaction_regulation <- reaction_regulation %>%
    transform_estimates(fun = log2, prefix = "log2_")

```

Significantly differentially regulated metabolic reactions are identified with function `identify_regulated()` 
called for our `actiData` object. 
The significant effects are defined by a numeric cutoff of raw or FDR-adjusted p-value and 
a numeric cutoff of fold-regulation estimate of reaction activity (identity scale!). 
The differential regulation status of the reaction, activated, inhibited or non-significant (ns), 
is stored in column `regulation` of the `reg` data frame in `actiData` objects.
With the code below, we'll define the significantly regulated metabolic reaction 
with FDR-adjusted p-value < 0.05 and at least 10% change (1.2-fold) in activity between 
the collagen subset of cancer samples:

```r

  reaction_regulation <- reaction_regulation %>%
      identify_regulated(p_type = "adjusted",
                         p_cutoff = 0.05,
                         fold_cutoff = 1.1)

```

Briefly inspected, our `actiData()` has now columns with $log_2$-transformed statistics, 
a columns storing the differential regulation status of metabolic reactions: 

```r

>   reaction_regulation
`actiData` object with activity estimates of 4742 in 94 subsystems

# A tibble: 4,742 × 17
   id      name  subsystem fold_reg  error lower_ci upper_ci p_value p_adjusted log2_fold_reg log2_error log2_lower_ci log2_upper_ci
   <chr>   <chr> <chr>        <dbl>  <dbl>    <dbl>    <dbl>   <dbl>      <dbl>         <dbl>      <dbl>         <dbl>         <dbl>
 1 R_13DA… 1,3-… beta-Ala…    1.46  0.130     1.26     1.80   0.0005   0.000797        0.547      -2.94          0.331        0.848 
 2 R_24_2… 24R-… Vitamin …    2.37  0.763     1.35     4.54   0.003    0.00438         1.24       -0.391         0.437        2.18  
 3 R_24_2… 24R-… Vitamin …    2.37  0.763     1.35     4.54   0.003    0.00438         1.24       -0.391         0.437        2.18  
 4 R_25VI… 1-al… Vitamin …    0.724 0.0583    0.617    0.850  0.0005   0.000797       -0.466      -4.10         -0.696       -0.235 
 5 R_25VI… 1-al… Vitamin …    0.724 0.0583    0.617    0.850  0.0005   0.000797       -0.466      -4.10         -0.696       -0.235 
 6 R_2AMA… 2-Am… Glycine,…    1.91  0.204     1.57     2.37   0.0005   0.000797        0.935      -2.29          0.651        1.25  
 7 R_2AMA… 2-Am… Methioni…    0.961 0.332     0.536    1.93   0.384    0.412          -0.0572     -1.59         -0.898        0.950 
 8 R_2HBO  2-Hy… Propanoa…    0.975 0.140     0.793    1.36   0.363    0.391          -0.0362     -2.83         -0.335        0.447 
 9 R_2OXO… 2-Ox… Lysine m…    0.916 0.0225    0.871    0.958  0.0005   0.000797       -0.126      -5.48         -0.199       -0.0614
10 R_34DH… 3,4-… Tyrosine…    2.19  0.520     1.52     3.90   0.0005   0.000797        1.13       -0.944         0.603        1.96  
# ℹ 4,732 more rows
# ℹ 4 more variables: p_type <chr>, p_cutoff <dbl>, fold_cutoff <dbl>, regulation <fct>
# ℹ Use `print(n = ...)` to see more rows

```

By counting of the regulation status in the `reg` data frame, 
we can see that roughly 1700 reactions were 
found to be significantly activated, while approximately 1800 reactions were inhibited 
in the collagen^high^ as compared with collagen^low^ subset of prostate cancer: 

```r 
>   reaction_regulation$reg %>% count(regulation)

# A tibble: 3 × 2
  regulation     n
  <fct>      <int>
1 activated   1711
2 inhibited   1831
3 ns          1200

```

Prior to any further analyses of differentially regulated reactions, we'll take a 
look at the very basic diagnostic of the model by plotting the prediction errors 
against predictions of fold-regulation of reaction activity. 
This is conveniently done with `plot()` method which returns a `ggplot` graphics 
with the requested scatter plot. 
The `fun = log2` argument specifies that we plot the estimates and their errors 
following a $log_2$ transformation. 

```r 

  est_error_plot <- reaction_regulation %>%
      plot(plot_type = "errors",
           fun = log2,
           plot_title = "Estimates and errors of reaction regulation",
           x_lab = "log2 fold-regulation estimate",
           y_lab = "log2 error")

```

In the plot, we immediately spot that errors inflate for large estimates of 
reaction activity - a point which, depending on your study question, may need 
discussion and changes of transcriptome quantification techniques. 
Our experience with `get_regulation()` function tells that the error inflation 
may be more efficiently addressed by modification of differential gene expression analysis 
methods than by tuning of arguments of `get_regulation()`.

---scatter plot of estimates and errors-----

#### Enrichment analyses for metabolic subsystems

Manual interpretation of roughly 3500 differentially regulated reactions 
is a daunting task; our experiences with general-purpose large language models/AI 
are rather disappointing. 
An "old-fashioned" but effective way to tackle such data is an enrichment analysis 
for metabolic subsystems within the sets of significantly activated and significantly 
inhibited reactions as compared with the whole reaction pool. 
To automate this step, we developed `suba()` function called for `actiData` objects. 
Please note that the function works only for objects with differential regulation status 
of the reactions - so make sure you to run `identify_regulated()` prior to calling 
`suba()`!
`suba()` offers two types of enrichment tests. 
By default, Fisher's exact test is performed (argument `type = "fisher"`). 
The alternative is random sampling of the entire reaction pool, and comparing the 
frequency of subsystems in these random samples with the frequency among the 
activated and inhibited reactions (`type = "random"`). 
With the code below, we'll conduct such random enrichment testing. 
`n_iter = 2000` argument specifies the number of random draws from the entire 
reaction pool:

```r

  subsystem_enrichment <- reaction_regulation %>%
      suba(type = "random",
           n_iter = 2000,
           .parallel = TRUE)

```

The result is a data frame with information on the tested reaction set 
(`reaction_set`: significantly activated or inhibited reactions), total number of 
reaction in the activated or inhibited reaction set (`n_total_reaction_set`), 
subsystem (`subsystem`), total number of reactions in the subsystem (`n_total_subsystem`), 
number of significantly activated or inhibited reaction in the subsystem (`n_intersect`). 
Magnitude of enrichment in the activated or inhibited reaction set as compared with 
the reaction pool is metered with odds ratio (OR, `or`) statistic with its confidence 
intervals (`lower_ci` and `upper_ci`).  
Finally the function returns raw p-values (`p_value`) and p-values corrected for 
multiple tests with the FDR method (`p_adjusted`): 

```r
>  subsystem_enrichment

# A tibble: 188 × 10
   reaction_set n_total_reaction_set subsystem              n_total_subsystem n_intersect    or lower_ci upper_ci p_value p_adjusted
   <fct>                       <dbl> <chr>                              <dbl>       <dbl> <dbl>    <dbl>    <dbl>   <dbl>      <dbl>
 1 activated                    1711 Alanine and aspartate…                11           5 1.39     0.857     6     0.355     0.495  
 2 activated                    1711 Alkaloid synthesis                     3           1 1.16     0.667     2     0.732     0.740  
 3 activated                    1711 Aminosugar metabolism                 25           7 0.857    0.571     1.6   0.268     0.427  
 4 activated                    1711 Androgen and estrogen…                31           6 0.607    0.438     1.17  0.036     0.094  
 5 activated                    1711 Arachidonic acid meta…                41          10 0.721    0.5       1.1   0.0725    0.162  
 6 activated                    1711 Arginine and Proline …                29          10 1.01     0.733     2.2   0.636     0.665  
 7 activated                    1711 beta-Alanine metaboli…                11           3 0.906    0.571     4     0.380     0.508  
 8 activated                    1711 Bile acid synthesis                   94          34 1.02     0.813     1.35  0.526     0.589  
 9 activated                    1711 Biotin metabolism                      8           0 0.297    0.167     1     0.026     0.0698 
10 activated                    1711 Blood group synthesis                 45          29 1.82     1.36      3     0.0005    0.00276
# ℹ 178 more rows
# ℹ Use `print(n = ...)` to see more rows

```

For our purposes, we'll define significant enrichment with FDR-adjusted p < 0.05 
and OR $\geq$ 2. 
This yields 16 significant subsystems: five subsystems significantly enriched 
among activated reactions and 11 subsystems significantly enriched within significantly 
inhibited reactions:

```r

 ## significantly enriched subsystems:
  ## OR >= 2 and pFDR < 0.05

  significant_subsystems <- subsystem_enrichment %>%
    filter(or >= 2, p_adjusted < 0.05)

```

```r

> significant_subsystems
# A tibble: 16 × 10
   reaction_set n_total_reaction_set subsystem              n_total_subsystem n_intersect    or lower_ci upper_ci p_value p_adjusted
   <fct>                       <dbl> <chr>                              <dbl>       <dbl> <dbl>    <dbl>    <dbl>   <dbl>      <dbl>
 1 activated                    1711 Chondroitin synthesis                 42          33  2.21     1.62     3.78  0.0005    0.00276
 2 activated                    1711 Glycosphingolipid met…                 7           6  2.38     1.75     7     0.013     0.0407 
 3 activated                    1711 Keratan sulfate synth…                59          58  2.73     2.11     4.21  0.0005    0.00276
 4 activated                    1711 Vitamin A metabolism                  29          21  2.04     1.47     4.4   0.0005    0.00276
 5 activated                    1711 Vitamin D metabolism                  13          10  2.15     1.22     5.5   0.0035    0.0150 
 6 inhibited                    1831 Biotin metabolism                      8           8  2.58     1.5      9     0.0015    0.00671
 7 inhibited                    1831 Cholesterol metabolism                47          37  2.05     1.52     3.17  0.0005    0.00261
 8 inhibited                    1831 Chondroitin sulfate d…                34          28  2.17     1.53     3.62  0.0005    0.00261
 9 inhibited                    1831 Citric acid cycle                     17          13  2.02     1.4      6.66  0.0015    0.00671
10 inhibited                    1831 Glutathione metabolism                15          13  2.27     1.4      4.67  0.0005    0.00261
11 inhibited                    1831 Keratan sulfate degra…                73          58  2.06     1.64     2.95  0.0005    0.00261
12 inhibited                    1831 N-glycan degradation                  14          11  2.08     1.5      6     0.003     0.0118 
13 inhibited                    1831 Oxidative phosphoryla…                 9           9  2.58     1.67    10     0.0005    0.00261
14 inhibited                    1831 Ubiquinone synthesis                   5           5  2.44     1.5      6     0.011     0.0369 
15 inhibited                    1831 Valine, leucine, and …                35          28  2.11     1.61     4.14  0.0005    0.00261
16 inhibited                    1831 Xenobiotics metabolism                20          19  2.48     1.54     5     0.0005    0.00261

```

We can plot the OR statistics and FDR-corrected p-values as a simple bubble plot 
with ad-hoc R code: 

```r

significant_subs_or_plot <- significant_subsystems %>%
    ggplot(aes(x = or,
               y = reorder(subsystem, or),
               size = -log10(p_adjusted),
               fill = reaction_set)) +
    facet_grid(reaction_set ~ .,
               scales = "free",
               space = "free") +
    geom_point(shape = 21) +
    scale_size_area(labels = function(x) signif(10^-(x), 2),
                    max_size = 4.5,
                    name = "pFDR") +
    scale_fill_manual(values = c("activated" = "firebrick",
                                 "inhibited" = "steelblue"),
                      name = "reaction set") +
    theme_classic() +
    theme(axis.title.y = element_blank()) +
    labs(title = "RECON metabolic subsystems, collagen hi vs low",
         x = "enrichment, OR")

```

----plot of OR and p-values-------

To explore these significantly enriched metabolic subsystems, we can tally 
their differentially regulated reactions by selecting the subsystems of interest 
in the `actiData` object and counting the differentially regulated reactions. 
To this end, we'll use functions `select()` and `count_regulated()`: 

```r

  reaction_regulation %>%
      select(subsystems = significant_subsystems$subsystem) %>%
      count_regulated

```
```r
# A tibble: 34 × 5
   regulation     n n_total percent subsystem                      
   <fct>      <int>   <int>   <dbl> <chr>                          
 1 activated    158     427   37.0  all                            
 2 inhibited    238     427   55.7  all                            
 3 activated      0       8    0    Biotin metabolism              
 4 inhibited      8       8  100    Biotin metabolism              
 5 activated      4      47    8.51 Cholesterol metabolism         
 6 inhibited     37      47   78.7  Cholesterol metabolism         
 7 activated      6      34   17.6  Chondroitin sulfate degradation
 8 inhibited     28      34   82.4  Chondroitin sulfate degradation
 9 activated     33      42   78.6  Chondroitin synthesis          
10 inhibited      3      42    7.14 Chondroitin synthesis          
# ℹ 24 more rows
# ℹ Use `print(n = ...)` to see more rows
```

Theses counts can be also displayed in a bar plot by calling `select()` and 
`plot()`. 
By setting the argument `type = "numbers"` of the `plot()` method, we specifically 
request a bar plot with counts (`scale = "counts"`) or percentages (`scale = "percent"`) 
of activated and inhibited reaction in the subsystems. 
Please note that total numbers of metabolic reactions in the subsystems are 
per default indicated in the Y axis: 

```r

  significant_subs_percent_plot <- reaction_regulation %>%
      select(subsystems = significant_subsystems$subsystem) %>%
      plot(plot_type = "numbers",
           scale = "percent",
           plot_title = "RECON metabolic subsystems, collagen hi vs low")

```

----bar plot of reaction numbers ------


### Calculation of activity scores of metabolic reactions with gene expression data


## References

[^1]: King ZA, Lu J, Dräger A, Miller P, Federowicz S, Lerman JA, Ebrahim A, Palsson BO, Lewis NE. 
BiGG Models: A platform for integrating, standardizing and sharing genome-scale models. 
Nucleic Acids Research (2016) 44:D515–D522. doi: 10.1093/NAR/GKV1049

[^2]: Gavai AK, Supandi F, Hettling H, Murrell P, Leunissen JAM, Van Beek JHGM. 
Using Bioconductor Package BiGGR for Metabolic Flux Estimation Based on Gene Expression Changes in Brain. 
PLOS ONE (2015) 10:e0119016. doi: 10.1371/JOURNAL.PONE.0119016

[^3]: Thiele I, Swainston N, Fleming RMT, Hoppe A, Sahoo S, Aurich MK, 
Haraldsdottir H, Mo ML, Rolfsson O, Stobbe MD, et al. 
A community-driven global reconstruction of human metabolism. 
Nat Biotechnol 2013 315 (2013) 31:419–425. doi:10.1038/nbt.2488

[^4]: Swainston N, Smallbone K, Hefzi H, Dobson PD, Brewer J, Hanscho M, Zielinski DC, 
Ang KS, Gardiner NJ, Gutierrez JM, et al. 
Recon 2.2: from reconstruction to model of human metabolism. 
Metabolomics (2016) 12. doi:10.1007/S11306-016-1051-4

[^5]: Heidegger I, Frantzi M, Salcher S, Tymoszuk P, Martowicz A, Gomez-Gomez E, 
Blanca A, Lendinez Cano G, Latosinska A, Mischak H, et al. 
Prediction of Clinically Significant Prostate Cancer by a Specific Collagen-related 
Transcriptome, Proteome, and Urinome Signature. 
Eur Urol Oncol (2024). doi:10.1016/J.EUO.2024.05.014
