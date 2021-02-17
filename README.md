# *noisyR*: Enhancing biological signal in sequencing datasets by characterising random technical noise #

High-throughput sequencing enables an unprecedented resolution in transcript quantification, at the cost of magnifying the impact of technical noise. The consistent reduction of random background noise to capture functionally meaningful biological signals is still challenging. Intrinsic sequencing variability introducing low-level expression variations can obscure patterns in downstream analyses.

The noisyR package comprises an end-to-end pipeline for quantifying and removing technical noise from HTS datasets. The three main pipeline steps are [i] similarity calculation across samples, [ii] noise quantification, and [iii] noise removal; each step can be finely tuned using hyperparameters; optimal, data-driven values for these parameters are also determined.

<img src="https://github.com/Core-Bioinformatics/noisyR/blob/master/docs/figures/workflow.png" width="400">

*Workflow diagram of the **noisyr** pipeline*

## Similarity calculation across samples ##

For the sample-similarity calculation, two approaches are available:

1. The **count matrix approach** uses the original, un-normalised count matrix, as provided after alignment and feature quantification; each sample is processed individually, only the relative expressions across samples are compared. Relying on the hypothesis that the majority of genes are not DE, most of the evaluations are expected to point towards a high similarity across samples. Choosing from a collection of >45 similarity metrics, users can select a measure to assess the localised consistency in expression across samples. A sliding window approach is used to compare the similarity of ranks or abundances for the selected features between samples. The window length is a hyperparameter, which can be user-defined or inferred from the data (supplementary methods 1). 
2. The **transcript approach** uses as input the alignment files derived from read-mappers (in BAM format). For each sample and each exon, the point-to-point similarity of expression across the transcript is calculated across samples in a pairwise all-versus-all comparison. 

The output formats for the two approaches are the same; the number of entries varies, since the count approach focuses on windows, whereas for the transcript approach we calculate a similarity measure for each transcript.

Main functions: *calculate_distance_matrices_counts()*, *calculate_distance_matrices_transcript()* 
Input preparation functions: *cast_matrix_to_double()*, *optimise_window_length()*, *cast_gtf_to_genes()* 

## Noise quantification ##

The noise quantification step uses the abundance-correlation (or other similarity measure) relation calculated in step i to determine the noise threshold, representing the abundance level below which the gene expression is considered noisy e.g. if a correlation threshold is used as input then the corresponding abundance from a (smoothed) abundance-correlation line plot is selected as the noise threshold for each sample. The shape of the distribution can vary across experiments; we provide functionality for different thresholds and recommend the choice of the one that results in the lowest variance in the noise thresholds across samples. Options for smoothing, or summarising the observations in a box plot and selecting the minimum abundance for which the interquartile range (or median) is consistently above the correlation threshold are also available. Depending on the number of observations, we recommend using the smoothing with the count matrix approach, and the boxplot representation with the transcript option.

<img src="https://github.com/Core-Bioinformatics/noisyR/blob/master/docs/figures/PCC_abn.png" width="500">

*Indicative plots of the Pearson correlation calculated on windows of increasing average abundance for the count matrix-based noise removal approach (left) and per exon for the transcript-based noise removal approach (right).*

Main function: *calculate_threshold_noise()* 
Visualisation function: *plot_distance_abundance()* 

## Noise removal ##

The third step uses the noise threshold calculated in step ii to remove noise from the count matrix (and/or BAM files). Genes/exons whose expression is below the noise thresholds for every sample are removed from the count matrix. The average noise threshold is calculated and added to every entry in the count matrix. This ensures that the fold-changes observed by downstream analyses are not biased by low expression, while still preserving the structure and relative expression levels in the data. If downstream analysis does not involve the count matrix, the thresholds obtained in step ii can be used to inform further processing and potential exclusion of some genes/exons from the analysis.

Main function: *remove_noise_method()* 

## Required packages ##

### CRAN (use *install.packages()*) ###
* utils
* grDevices
* tibble
* dplyr
* magrittr
* ggplot2
* philentropy
* doParallel
* foreach

### Bioconductor (use *BiocManager::install()*) ###
* preprocessCore
* IRanges
* GenomicRanges
* Rsamtools
