---
output: 
  md_document:
    variant: markdown_github
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

# Biodiv

<h5 align="right">

Latest version: 2023-08-18

</h5>

<font color="394CAE">

<h3 color="394CAE" style="font-weight: bold">

Introduction to Biodiv (R package): Excerpt from Biodiv User's Guide

</h3>

</font> <br>

<h5><b> Wesley Neves</b> <br><br></h5>

<br> Biodiv (Package for calculating Taxonomic Diversity, Functional Diversity, and Functional Redundancy) is an R package, avaliable in [Github](https://github.com/wesneves/Biodiv.git). In this document, here provide a quick introduction demonstrating how to run the package `Biodiv` (Taxonomic Diversity, Functional Diversity and Functional Redundance). `Biodiv` has several main functions: `FDsingle`, `FDchao`, `Div` and `PlotDiv`.

The `Biodiv` is an package the extension of function [FunD] ([https://github.com/AnneChao/FunD.git)](https://github.com/AnneChao/FunD.git)) (Chao et al , 2019) `Biodiv` focuses on three measures of Hill numbers of order q: species richness (`q = 0`), Shannon diversity (`q = 1`, the exponential of Shannon entropy) and Simpson diversity (`q = 2`, the inverse of Simpson concentration) and extend Hill numbers to three dimensions under Hill-Chao family frame work (Chao et al., 2019): Taxonomic diversity (TD), Functional diversity (FD) and Functional redundance (FR).

The functional diversity of order q at level $\tau>0$ as the attribute diversity of order q:

$${}^q F D(\Delta(\tau))=(\sum_{i=1}^8 v_i(\tau)(\frac{a_i(\tau)}{n_{+}})^q)^{1 /(1-q)}\\ =(\sum_{i=1}^S p_i(\sum_{j=1}^S[1-f(d_{i j}(\tau))] p_j)^{q-1})^{1 /(1-q)}$$

The Level of threshold distinctiveness is:

$$
\begin{equation} \Delta(\tau) \equiv\left[d_{i j}(\tau)\right]=\left[\min \left(d_{i j}, \tau\right)\right] . \end{equation}
$$

The value tau is obtained through the index Rao's quadratic entropy Q, which represents the mean functional distance between any two individuals randomly selected from the assemblage (hence accounting for species abundance).

The diversity functional is obtained through the tau value minimun and medium according to RAO Q index, and the functional redundancy was based on a methodological modification of the work by (Ricotta et al 2016), in which functional redundancy was calculated as the relative difference between the Rao's Q index ratio and the Simpson index.

$$FR_{R}=1-\frac{Q(X)}{D(X)}=1-\frac{(\sum_{j=1}^S p_j \sum_{i \neq j}^S(d_{i j} p_i))}{(\sum_{j=1}^S p_i \sum_{i \neq j}^S p_i)}$$

The result is expressed within the interval [0, 1]. However, since these analyses do not strictly adhere to several desirable mathematical properties as described in (Jost 2009), several authors (Chao et al, 2014, 2019; Jost, 2006; Leinster & Cobbold, 2012) demonstrated that these inconveniences are overcome by converting the raw entropy values into what is known as effective species. Yet, these analyses solely focus on the concept of taxonomic diversity. Therefore, authors (Leinster & Cobbold, 2012 and Chao et al 2019) suggested that Hill numbers should take into account all equivalent forms of species diversity indices from the perspectives of functional and phylogenetic diversity as well. Both authors developed formulas known as generalized Hill numbers (Chao et al, 2019, Leinster & Cobbold, 2012). For instance, Leinster and Cobbold relate various dissimilarity-based indices to one another, based on the weight they assign to common functional groups. Chao and collaborators, on the other hand, not only considers dissimilarities between pairs of species but also weighs these values against a threshold of pairwise species distinction (Tau), thereby making the analysis more meaningful in representing the community from a functional perspective. Since in Chao's work the authors argue that the results of the methodology proposed by Leinster & Cobbold will always lead to underestimated measurements of biological diversity.

The function `FDchao` calculate Functional Diversity of N sites for various values of tau and q.

The function `Div` calculate Functional Redundance blending the approach from the work of Ricotta et al 2016 with that of Chao et al 2019, the functional redundancy obtained through this new approach is essentially the relative difference between the ratio of functional diversity weighted by the mean value of Tau and the biological diversity weighted by the minimum value of tau (corresponding to classic taxonomic diversity, which can be derived using conventional Hill numbers).

$$
FR=1-\frac{{}^qFD(\Delta(\tau_{med}))}{{}^qFD(\Delta(\tau_{min}))}=1-\frac{\left(\sum_{i=1}^S v_i(\tau)\left(\frac{a_i(\tau)}{n_{+}}\right)^q\right)^{1 /(1-q)}}{\left(\sum_{i=1}^{S} p_{i}^{q}\right)^{1/(1-q)}}
$$

The value of the generalized index for functional diversity, when weighted by the minimum value of the functional distance matrix, results in taxonomic diversity since all species are grouped into individual clusters, as no pair of species has a distance smaller than the smallest value in the distance matrix. Therefore, no pairwise species value is smaller than the minimum value of tau.

The result is saved in a data.frame object necessary for plotting graphs of taxonomic, functional, and redundancy diversity.

## Install Biodiv in R

```{r}

## install the latest version from github
install.packages('devtools')
library(devtools)
install_github('wesneves/Biodiv')
```

### MAIN FUNCTION: FDchao(), Div(), PlotDiv()

``` r
FDchao(data, distance, tau, q, boot)
```

This function calculate Functional Diversity of N sites for various values of tau and q

``` r
Div(data)
```

Data.frame for value of diversity required to plot the results. This function prepares a "data.frame" object where it takes 5 variables, namely (qEix, TauMin, TauMed, RedFunChao, Color). These variables are required to plot the graphs (using the PlotDiv function) depicting taxonomic, functional, and functional redundancy diversity results.

``` r
plotDiv(data, tog, cap)
```

Function to plot the value of diversity. This is a function to plot the results of the FDchao function filtered by the Div function, generating three graphs: Taxonomic Diversity, Functional Diversity, and Functional Redundancy.

The parameter `data` must necessarily be an object containing the outcome of the `Div` function to prepare the graph, `tog` is an abbreviation of the word "together" and is a logical object. If it is set to `TRUE`, the function will plot the three graphs together (side by side), `cap` is an abbreviation of the word "capitions" and is a logical object if it is set to `TRUE`, the function will ask you to set the labels for the graph.

## Data Format:

The type of data are supported is Individual-based abundance data Input data for each assemblage/site include samples species abundances in an empirical sample of n individuals ("reference sample"). When there are N assemblages, input data consist of an S by N abundance matrix, or N lists of species abundances.

## Example:

```{r}
abunlist <- as.list(abundm) #list abundance for species
```

```{r}
head(abundm, 3)
                      pv        ia        tm ga        cq
A_brasiliensis 0.6666667 0.0000000 1.3333333  1 0.6666667
Bledius        0.0000000 0.3333333 0.0000000  0 0.3333333
Coelopa        0.0000000 0.0000000 0.6666667  0 0.0000000
```

```{r}
ktab_cat <- ade4::ktab.list.df(list(trait)) #Trait is an matrix for categorical data of macrobenthic species 
dist_cat <- ade4::dist.ktab(ktab_cat, type = "N")
distmat <- as.matrix(dist_cat)

library(Biodiv)
library(ggplot2)
library(magrittr)

FDresul <- FDchao(abunlist, distmat, seq(min(distmat[distmat>0]), max(distmat), length.out = 25), seq(from = 0, to = 2, length.out = 25), 50)
```

```{r}
head(FDresul$forq, 3)
           q estimate   s.e.    LCL     UCL  tau site
1 0.00000000  11.0000 1.3157 8.4213 13.5787 dmin   pv
2 0.08333333   9.9450 1.0299 7.9264 11.9636 dmin   pv
3 0.16666667   8.9997 1.0606 6.9209 11.0785 dmin   pv
```

The variable `estimate` represents the outcomes of diversity indices applied with a functional diversity weight (tau value).

```{r}
FDfilt <- Div(FDresul)
```

```{r}
head(FDfilt, 3)
        qEix  TauMin TauMed RedFunChao Color
1 0.00000000 11.0000 7.5219  0.3161909    pv
2 0.08333333  9.9450 6.9605  0.3001006    pv
3 0.16666667  8.9997 6.4687  0.2812316    pv
```

We used the following command to obtain estimates of the minimum and average value of the parameter Tau.:

Minimum Value:

```{r}
min(distmat[distmat>0]) #min value
```

Medium Value:

```{r}
abunlist %>%
    lapply(FUN = function(x) x/sum(x)) %>%
    do.call(cbind, .) %>%
    apply(1, mean) %>%
    {sum((.%*%t(.))*distmat)} #med value
```

Plotting the graphs for Functional Diversity, Taxonomic Diversity, and Functional Redundancy:

Example 1: Plotting the three graphs side by side with predefined axis labels, title of the graphs:

```{r}
abundm <- read.csv("~/Documentos/Trabalhos R/Dissertation_UFPE/Dados/abundm.csv", row.names=1)
abunlist <- as.list(abundm)

trait <- read.csv("~/Documentos/Trabalhos R/Dissertation_UFPE/Dados/trait21.csv", row.names=1)
ktab_cat <- ade4::ktab.list.df(list(trait))
dist_cat <- ade4::dist.ktab(ktab_cat, type = "N")
distmat <- as.matrix(dist_cat)
library(Biodiv)
library(ggplot2)
FDresul <- FDchao(abunlist, distmat, seq(min(distmat[distmat>0]), max(distmat), length.out = 25), seq(from = 0, to = 2, length.out = 25), 50)
FDfilt <- Div(FDresul)
plotDiv(FDfilt, tog = T, cap = F)
```

```{r}
# Example 2: Plotting the three graphs side by side while modifying the axis labels and the title of the graphs.
plotDiv(FDfilt, tog = T, cap = T)

# Example 3: Plotting one graph with predefined axis labels, title of the graph
plotDiv(FDfilt, tog = F, cap = F)

# Example 4: Plotting one graph while modifyinh the axis labels and the title of the graphs.
plotDiv(FDfilt, tog = F, cap = T)
```
