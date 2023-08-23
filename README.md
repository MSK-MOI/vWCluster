# vWCluster
vWCluster is a network-based clustering method based on the vector-valued Wasserstein distance derived from optimal mass transport (OMT) theory as presented in [vWCluster: Vector-valued optimal transport for network based clustering using multi-omics data in breast cancer](https://doi.org/10.1371/journal.pone.0265150)

vWCluster treat mRNA expression, DNA copy number alteration, and DNA methylation data a multi-layer super graph structure and a novel extension of optimal tranport, vector-valued OMT is used. 

A sample of this data for 1904 breast cancer (METABRIC) patients and 3417 genes is given in the folder "data". 

## Requirements:
1. MATLAB
2. Install CVX: http://cvxr.com/cvx/doc/install.html prior to running this function.

## In this repository:

1. main finds the pairwise Wasserstein distances between samples. 

2. Invarinat_3: creates the integative measure of (mRNA, CNA and Methylation) for a given sample.

3. dist_cvx: individual vector-valued OMT distance calcualtion

4. Hier_clustering build a hierachical clustering of the samples using the Wasserstein distances and then run Kaplan-Meier analysis

5. MatSurv: Kaplan-Meier analysis, download from https://www.mathworks.com/matlabcentral/fileexchange/64582-matsurv
