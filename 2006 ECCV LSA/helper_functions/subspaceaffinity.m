function d=subspaceaffinity(base1, base2)
%function d=subspacedistance(base1, base2)
%Compute an affinity measure between two subspaces represented by the ORTHONORMAL
%basis BASE1 and BASE2
theta=subspaceangle(base1,base2);
d=exp(-sum((sin(theta).^2)));
