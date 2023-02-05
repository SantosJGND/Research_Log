# Stats Lab

*Exploring and classifying genetic data.*

In recent years the size and quality of biological data sets have increased at a tremendous rate. It is fair to say that we 
have not yet reached the plateau of that development, so that students and researchers alike will be facing even larger data 
sets in the near future.

As the amount of information increases, so must the methods we use to explore, analyse and even present that data. 
The field of population genetics is undergoing a profound paradigm shift, as researchers move from exploiting the sparse, 
hard-won data sets of old to finding ways to understand and encompass the swaths of data their laboratories turn in. 

In recent years i have had the privilege of working on one of the largest genetic data sets to date (2018). In summary, this 
project involved the description of genetic structure across thousands of possible sites. I would like to leave 
here some of the more useful insight and tools that resulted. This repository consists in a series of Jupyter notebooks, 
each exploring some method or applications i found particularly useful.

Notebooks are organised in order of increasing complexity and, in most cases, of chronological development. This development was driven in 
greater part by the need to study the relation between statistics that describe structure in data. Metrics to describe structure vary
across fields of research. Here, the focus is on the description of genetic data and on phased SNP data in particular. 

Use NBviewer directly to explore the [notebook library](https://nbviewer.jupyter.org/github/SantosJGND/Genetic-data-analysis/tree/master/). 
[NBviewer](https://nbviewer.jupyter.org/) supports the 2D and 3D images produced in the scripts (copy-paste the url of a notebook onto the tab provided).


## I. Binary samples. Population distribution and identification.

- Notebooks: 1 - 5

Structure in data is either used to study real-world processes influencing the distribution of observed variables or as a basis for
prediction. From a mathematical point of view these two goals are indistinguishable. In either case the focus of scientific research is
on the description of the laws that govern the generation of new data. This description comes in the form of the combined distributions
of all variables available to describe the subject of study. 

In the study of the variation of binary data, the simplest assumption is that of a binomial probability of observation at the smallest 
unit of measurement  available. In the case of genetic data, this unit is the single nucleotide polymorphism (SNP).

Througout the notebooks provided, individuals are simulated as samples of *L* binary markers (0, 1, .., i), coded 0 and 1, of frequency *pi* 
within a given population *K*. Populations are simulated as multivariate Bernoulli variables *Sk*, where each variable *Ski* represents the 
binomial probability of an event. Allele frequencies *Sk* are drawn from the Beta distribution (Tataru *et al.* 2017, 
see also Jiang and Cockerham 1990, Sawyer and Hartl 1992 and Williamson *et al.* 2004). Population-specific allele frequency vectors 
of size *L* are sampled from the Beta distribution at various combinations of mean and variance of this distribution (**1**)
To the extent that individual variables are assumed to bear the same weight, principal component analysis (PCA) presents an intuitive
summarisation of population samples (see notebook **8** in the next section). In this context, given a sufficient number of samples, a 
basic description of a population entity is the probability density function of its samples in PCA feature space. 

In Notebooks **2**, **3**, **4** and **5**, the use of kernel density estimates ([KDE](https://scikit-learn.org/stable/modules/density.html)) 
for the characterisation and assignment of individual samples is explored. 

- [1. Generating samples](https://nbviewer.jupyter.org/github/SantosJGND/Genetic-data-analysis/blob/master/Notebooks/1.%20Generating_haplotypes.ipynb).
- [2. Local classification](https://nbviewer.jupyter.org/github/SantosJGND/Genetic-data-analysis/blob/master/Notebooks/2.%20Local_classification.ipynb).
- [3. Mislabelling](https://nbviewer.jupyter.org/github/SantosJGND/Genetic-data-analysis/blob/master/Notebooks/3.%20Mislabelling.ipynb).
- [4. Outlier material](https://nbviewer.jupyter.org/github/SantosJGND/Genetic-data-analysis/blob/master/Notebooks/4.%20X-material.ipynb)
- [5. Visualizing KDE](https://nbviewer.jupyter.org/github/SantosJGND/Genetic-data-analysis/blob/master/Notebooks/5.%20Visualizing%20KDE.ipynb)

## II. Population structure.

Notebooks: 6 - 11

The study of the relation of different descriptors to population structure requires the analysis of samples from as many different populations
as possible. In this section we begin by dealing with organising simulations so as to control for correlation (notebook **6**). We then build on
the infrastructure created to explore the impact of population distance on classification accuracy (same notebook) and introduce **intermediate 
classes** as a way to circumvent the problem of increased assignment error at lower distances ([Fig. 1](https://github.com/SantosJGND/Stats_Lab/blob/master/Complementary_data/Ideo_order_Fst_example.png)). In notebook **7** we introduce
varying population sizes and numbers, and explore their impact on the projection of samples in PCA space. Notebook **8** extends this study and provides 
a real life example.

Notebooks **9a-b** and **10** explore a particular way of studying population genetic structure. The focus of these notebooks is on capturing the
degree of overlap between distributions. The description explored, *overlap*, is not one of distance or variance (like *Fst* or ANOVA), but of the extent
of distribution space significant to more than one population. In terms of classification, this amounts to a study of error.

**The importance of the overlap measure**

In notebooks **6** through **8** we explored the relation between correlation, distance in PCA space and classification accuracy. We learned that for 
two populations this relation is straighforward: error rate presents a gaussian decline with increasing distance (*Fst*) but is positively and linearly
related to the overlap measure. The reason why this second finding is important is that natural populations are not always neat little blobs of gaussian
variation. The assumption of this ideal scenario is in fact one limitation in existing approaches in dealing with complex data sets. 

Natural populations can present an important degree of internal substructure. In cases where internal variation does not exceed the variance between 
reference clusters estimates of deviation from internal centroids can still be informative. However, if substructure is the result of exogenous introgression, 
then the accompanying increase in internal variance can complicate the assignment of individual samples. Enter KDE and the study of overlap. Through KDE
analyses become limited to those clusters of variation present in a given data set, all isolated samples being assigned as outliers 
(see notebook **4**, section I). It is possible for samples from the same reference populations to be unevenly distributed among local 
clusters. It is further possible for local clusters to be populated by samples from different populations. In this context, we need to know how 
the relative presence of reference samples at a given cluster will affect their assignment through KDE. However, because we are now dealing with structured populations
we can no longer rely on the *Fst* measure to interpret variation in error rate. This is the reason for the development of the overlap measure.

In notebook **9b** we explore the impact of asymmetric overlap of structured populations on both pure and intermediate classification (see notebook **6**).

- [6. A link to Fst](https://nbviewer.jupyter.org/github/SantosJGND/Genetic-data-analysis/blob/master/Notebooks/6.%20A%20link%20to%20FSTs.ipynb)
describes the use of PCA to organise allele frequency vectors by genetic proximity, as measured in *Fst* 
([Fig. 2](https://github.com/SantosJGND/Stats_Lab/blob/master/Complementary_data/Ideo_comparison.png)) and introduces intermediate classes.

- [7. Controlling for Size](https://nbviewer.jupyter.org/github/SantosJGND/Stats_Lab/blob/master/Notebooks/7.%20Controlling%20for%20size.ipynb) 
describes the relation between the proximity of source populations and euclidian distances of sample projections in PCA space.
The impact of sample size is studied.

- [8. Machine learning for Genomics](https://nbviewer.jupyter.org/github/SantosJGND/Stats_Lab/blob/master/Notebooks/8.%20Machine%20Learning%20for%20Genomics.ipynb)
describes the use of PCA to infer genetic distances measured in *Fst*.

- [9. Conditional variation](https://nbviewer.jupyter.org/github/SantosJGND/Genetic-data-analysis/blob/master/Notebooks/9.%20Conditional_variation.ipynb) 
and [9b. Conditional Var, but with Pops](https://nbviewer.jupyter.org/github/SantosJGND/Genetic-data-analysis/blob/master/Notebooks/9b.%20Conditional_Var%2C%20but%20using%20Pops.ipynb) 
describe the use of KDE to extract a measure of distribution overlap. 

- [10. Indexing MS](https://nbviewer.jupyter.org/github/SantosJGND/Genetic-data-analysis/blob/master/Notebooks/10.%20Indexing_MS.ipynb) 
and [11. Visualising overlap](https://nbviewer.jupyter.org/github/SantosJGND/Genetic-data-analysis/Notebooks/blob/master/11.%20Visualising%20Overlaps.ipynb) 
explore the behaviour of Mean Shift classification under different patterns of structure.


## III. Admixture.

Notebooks 12 - 13

This section extends on the analyses in Notebooks 2 through 5. The use of the KDE of known populations for the characterisation 
of samples in scenarios containing outliers and cases of mislabelling is explored. However, here we will build on the infrastructure
developed in Section II to compare classification output across structure and admixture scenarios.

- [12a. Admixture simulations](https://nbviewer.jupyter.org/github/SantosJGND/Genetic-data-analysis/blob/master/Notebooks/12a.%20Admixture%20simulations.ipynb) 
introduces individual IDs and admixture proportions to combine admixture and genetic structure ([Fig. 3](https://github.com/SantosJGND/Stats_Lab/blob/master/Complementary_data/3way_example.png)).

- [12b. Infering Admixed](https://nbviewer.jupyter.org/github/SantosJGND/Stats_Lab/blob/master/Notebooks/12b.%20Inference_Admixed.ipynb) 
KDE-based assignment of admixed scenarios.

- [13. Outliers and genetic structure](https://nbviewer.jupyter.org/github/SantosJGND/Stats_Lab/blob/master/Notebooks/13.%20Outliers_genetic_structure.ipynb) 
focuses on the identification of samples from non-reference sources under different structure parameters ([Fig. 4](https://github.com/SantosJGND/Stats_Lab/blob/master/Complementary_data/Supplemental_Figure_S11.png)).


## IV. Naked structure.

Throughout the previous sections the exploration of genetic structure was made through the study of local structured correlations. The focus was on the 
characterisation of samples relative to the distributions of supervised and unsupervised groups. However, while the use of KDE and of Mean Shift through it 
provides an accurate description these distributions, and allows for the use of outliers, it does so at the expense of one important element in the 
description of population genetic data: distance.

Except in the ideal scenario of normal unadmixed populations, the *p*-values derived from KDE cannot be used to infer the actual distances that separate 
reference groups. Even in that ideal case, our use of kernels would make the inference of distances beyond the threshold of 0.03 *Fst*s extremely imprecise. 
Why was this limitation allowed?

In the study of natural populations, the likelihood of polyphyly and recurrent genetic exchange, with possible confounding effects on population 
genetics estimates, is high. This is what led us, in section II, to rely on a measure of *overlap* to study error rate in the classification of samples
from structured populations. In the context of the analysis of genomic data, where most analyses are automated, the problem presented by these occurrences for 
the calculation of genetic distances is that we might be wrong in assuming the identity of the elements we are calculating distances between. This can introduce 
a considerable amount of variation in our measurements. More importantly, this factor can seriously impact any classification of individual 
samples reliant on distance measures. This led to the decision of analysing local correlation and distance measurements separately.

More than the previous sections, this section is about exploratory analysis. We will begin by exploring methods to extract and visualise distributions of distances
between populations. Once again relying heavily on plotly tools, notebook **14** explores how different patterns of genomic structure impact the
distribution of distances of one focal population to the remainder. notebook **15** shows how this analysis can be guided by the supervised classificaiton
of samples using KDE. In notebook **16** we develop the targeted analysis of distance to one of structure, treating multiple populations together.


- [14. Target distances](https://nbviewer.jupyter.org/github/SantosJGND/Stats_Lab/blob/master/Notebooks/14._Target_distances.ipynb). Study of differentiation between a target
and a control group. Introduces updated software for the manipulation of structure across data sets.

- [16. Structure vertices](https://nbviewer.jupyter.org/github/SantosJGND/Stats_Lab/blob/master/Notebooks/16._Structure_vertices.ipynb). Study the combination of distances
that comprise a population structure. Explore use for prediction.


## Folders and other directories.

- **/Notebooks** jupyter notebooks.

- **/Dmatrices** and **/Complementary_data** hold Darwin format files for notebook **8** and figures linked to in readme respectively.

- [/LAI_interface](https://github.com/SantosJGND/Stats_Lab/tree/master/LAI_interface): analysis of simulated samples using Local Ancestry Inference software 
([Fig. 5](https://github.com/SantosJGND/Stats_Lab/tree/master/LAI_interface/Supplemental_Figure_S6.png)).

- [Digits](https://github.com/SantosJGND/Digits): an aside on controlling for admixture proportions for many individuals as a function of 
proxy genomic position ([Fig. 6](https://github.com/SantosJGND/Digits/blob/master/Ideo_step__OutlierTh0.0001_Z1.png)).


## References

- Jiang CJ and Cockerham CC. 1987. Use of the multinomial Dirichlet model for analysis of subdivided genetic populations. Genetics **115**: 363-366.

- Sawyer SA and Hartl DL. 1992. Population genetics of polymorphism and divergence. Genetics **132**: 1161-1176.

- Tataru P, Simonsen M, Bataillon T, Hobolth A. 2017. Statistical inference in the Wright-fisher model using allele frequency data. *Syst. Biol.* **66**:e30–e46.

- Williamson S, Fledel-Alon A and Bustamante CD. 2004. Population genetics of polymorphism and divergence for diploid selection models with arbitrary dominance. Genetics **168**: 463-475.
