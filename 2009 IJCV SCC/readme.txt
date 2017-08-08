This toolbox contains the Matlab implementation of the Spectral Curvature Clustering algorithm 
(by Guangliang Chen and Gilad Lerman, IJCV 2009). 

The main files are scc.m and lscc.m (the rest are supporting files). 
The former (scc) can be applied to segment general affine subspaces (including linear subspaces), 
while the latter (lscc) is restricted to linear subspaces. 

The simplest way to use them is as follows:

[sampleLabels, averageL2Error] = scc(X,d,K);
[sampleLabels, averageL2Error] = lscc(X,d,K).

where,

X: N-by-D data matrix 
d: intrinsic dimension
K: number of clusters

sampleLabels: cluster labels of the data points
averageL2Error: averaged L2 error of the detected clusters.

For more detailed description, as well as choices of optional parameters, 
please type in the Matlab command window

help scc

If you have any questions please email Guangliang Chen at glchen@math.duke.edu 
or Gilad Lerman at lerman@umn.edu.