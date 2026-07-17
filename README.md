# COTRA
This repository contains code and tests for the COTRA pipeline (CRISPR On-Target Ranking Architecture), a framework for identifying the single on-target instance in single genome-editing experiment data and for assessing dataset difficulty in distinguishing on-target from off-target instances.
# Dataset description
FINAL_RESULTS_1732_4894_MMs.tsv contains the counted edits done on Hek293T cell line with the help of two editing agents: spCas9 nuclease and BEs and two gRNA (+PAMs): TGCTCGAGTGGGTCCCCGTGAGG and GAGGACGAGATGTAAGAGGCTGG targeting the 1732 and 4894 transcript position of the COL7A1 gene, respectively. 
# Description of the particular columns
*chrom:* the contig which the edit is located on,  
*start:* the coordinate of the beginning of the edit,  
*end:* the coordinate of the end of the edit,  
*strand:* orientation of the edited strand (+/-),  
*1732_1.breakends.noEnds.FILTERED : Neg_ctrl_1h_DNA1.breakends.noEnds.FILTERED:* number of edits in particular samples, calculated according to the INDUCE-seq protocol,  
*TOTAL_SUMS:* the rowwise sum of edits in a given loci. One locus is either one edited nucleotide or a group of nucleotides with pairwise distances smaller than 5nt apart,  
*cluster_id:* unique ID of the edited locus (see description of the *TOTAL_SUMS* column,  
*clustered_sequence:* nt sequence of the edited locus (either one or more nt long),  
*centered_sequence:* flanking nt sequence of the edited locus, 100nt long on both sides of the edit,  
*strand_collapsed_cluster_id:* unique ID of the edited locus, IRRESPECTIVELY of the strand orientation (i.e. same coordinates on opposite strands clustered together),  
*strand_collapsed_cluster_total:* the rowwise sum of edits of a given inter-strand cluster,  
*hit_orientation_1732:* the strand orientation of the optimal alignment of the 1732 gRNA+PAM sequence to the flanking sequence of the edit in the sense of optimizing local sliding-window Hamming distance,  
*best_start_1732:* the starting coordinate of the optimal alignment of the 1732 gRNA+PAM sequence, according to criterion described above,  
*best_end_1732:*  the ending coordinate of the optimal alignment of the 1732 gRNA+PAM sequence, according to criterion described above,  
*best_substring_1732:* the subsequence of the edit flanking sequence, with the length equal to the length of gRNA+PAM sequence, giving the lowest edit distance to gRNA as measured by Hamming distance,  
*mismatch_count_1732:* edit distance between optimal flanking subsequence (see description above) and gRNA+PAM sequences, as measured by Hamming distance,  
*hit_orientation_4894 : mismatch_count_4894:* same as *hit_orientation_1732 : mismatch_count_1732* but for 4894 gRNA+PAM;


# Modules description
**00_setup.R**: configuration of the working environment. Installing necessary packages, if needed;  
**01_io.R**: reading in the input data in a tabular form;  
**02_sampling.R**: splitting the input data into particular samples;  
**03_features.R**: extracts from the edit flanking sequence the features used subsequently to train the model. Among them:
 - GC content,
 - AT/GC skews,
 - CpG count,
 - k-mer frequencies (from order=1 to order=4);

All the itemized features gives a 344-dimensional feature space. This feature space, jointly with the edit occurrences per sample forms a **feature matrix**;  
**03b_preprocessing.R**: performs early preprocessing of the **feature matrix**. This preprocessing is aggregated into 3 levels:
 - *none*: no preprocessing is performed,
 - *basic*: includes: zero and near-zero variance removal, as well as removal of exact linear combinations. The last step is especially of importance because of compositional nature of k-mer representation, which induces at least one column being an exact linear combinations of others per k-mer order;
 - *aggressive*: includes all preprocessing steps from the *basic* level as well as a removal of highly correlated features (>0.95 as measured by the Pearson correlation);

**04_routing_blocks_hpca2.R**: This module constructs blocks of variables for hPCA-based dimensionality reduction.

First, it computes a set of statistics describing the sparsity and distribution of features in the feature matrix, including the Gini coefficient, feature entropy, fraction of zeros per feature, support, variance, standard deviation, and mean.

Based on these statistics, feature-specific cutoffs are defined to determine how similarity between features should be measured. Depending on feature density, similarity is computed using one of three metrics: Spearman correlation for dense features, the Hellinger coefficient for moderately dense features, and the Jaccard metric for sparse features.

As a result, features are effectively partitioned into three separate spaces, within which hierarchical clustering is performed independently to divide features into blocks.

Finally, clusters smaller than a user-defined threshold are iteratively reassigned to the most similar larger clusters until all clusters satisfy the minimum size requirement. Cluster similarity is assessed using the RV coefficient.

# Dependencies
"Biostrings", "ggplot2", "FactoMineR", "dynamicTreeCut", "multiblock", "e1071", "isotree", "ineq", "entropy", "vegan", "mclust", "NMF", "caret", "RGCCA", "Rdimtools", "fastICA"
