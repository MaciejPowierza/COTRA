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



# Modules description
TD
