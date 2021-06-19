# vWCluster
vWCluster is a network-based clustering method based on the vector-valued Wasserstein distance derived from optimal mass transport (OMT) theory as presented in https://www.biorxiv.org/content/10.1101/2021.06.17.448878v1.

vWCluster treat mRNA expression, DNA copy number alteration, and DNA methylation data a multi-layer super graph structure and a novel extension of optimal tranport, vector-valued OMT is used. 

A sample of these data for 100 breast cancer (tcga) patients and 290 genes (intersection of HPRD and OncoKB ) is given in the file "Data".

Requirements:
1. MATLAB
2. Install CVX: http://cvxr.com/cvx/doc/install.html prior to running this function.

In this repository:

1. main finds the pairwise Wasserstein distances between samples. 
```
Invarinat_3: creates the integative measure of (mRNA, CNA and Methylation) for a given sample.
```
dist_cvx: individual vector-valued OMT distance calcualtion

2.  Hier_clustering build a hierachical clustering of the samples using the Wasserstein distances and then run Kaplan-Meier analysis
```
MatSurv: Kaplan-Meier analysis, download from https://www.mathworks.com/matlabcentral/fileexchange/64582-matsurv
