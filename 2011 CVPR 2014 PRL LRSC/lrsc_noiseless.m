%--------------------------------------------------------------------------
% C = lrsc(A,tau)
% Low Rank Subspace Clustering algorithm for noiseless data lying in a 
% union of subspaces
%
% C = argmin |C|_* + tau/2 * |A - AC|_F^2 s.t. C = C'
%
% A: clean data matrix whose columns are points in a union of subspaces
% tau: scalar parameter 
%--------------------------------------------------------------------------
% Copyright @ Rene Vidal, November 2012
%--------------------------------------------------------------------------

function C = lrsc_noiseless(A,tau)

if nargin < 2
    tau = 100/norm(A)^2;
end

[~,S,V] = svd(A,0);
lambda = diag(S);
r = max(sum(lambda > 1/sqrt(tau)),1);
C = V(:,1:r) * (eye(r) - diag(1./(lambda(1:r).^2)/tau)) * V(:,1:r)';

