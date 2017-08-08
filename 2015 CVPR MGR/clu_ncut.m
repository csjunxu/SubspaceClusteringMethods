function [idx] = clu_ncut(L,K)
% this routine groups the data X into K subspaces by NCut
% inputs:
%       L -- an N*N affinity matrix, N is the number of data points
%       K -- the number of subpaces (i.e., clusters)
L = abs(L)+abs(L');

D = diag(sum(L,2).^(-1./2));
L = eye(size(L,1)) - D*L*D;
[U,S,V] = svd(L);

V = U(:,end-K+1:end);
idx = kmeans(V,K,'emptyaction','singleton','replicates',10,'display','off');
idx = idx';