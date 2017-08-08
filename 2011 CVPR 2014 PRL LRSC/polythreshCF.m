%--------------------------------------------------------------------------
% [A,C] = polythreshCF(D,tau,alpha)
% Low Rank Subspace Clustering algorithm for data lying in a union of 
% subspaces and contaminated with outliers
%
% min |C|_* + tau/2*|A-AC|_F^2 + alpha/2*|Delta-A|_F^2  s.t. C = C'
%
% C = affinity matrix
% A = clean data matrix whose columns are points in a union of subspaces
% Delta = data matrix whose columns are points in a union of subspaces
%         (without outliers)
% tau = scalar parameter 
% alpha = scalar parameter 
%--------------------------------------------------------------------------
% Copyright @ Rene' Vidal
% Edited by Paolo Favaro, January 2013
%--------------------------------------------------------------------------
function [A,C] = polythreshCF(D,tau,alpha)
% compute polynomial thresholding operator by using
% an approximate closed form formula

[Dx,Dy] = size(D);
[U,S,V] = svd(D);

if (3*tau < alpha)
    sigma = (alpha + tau)/alpha/sqrt(tau);
else
    sigma = sqrt((alpha + tau)/alpha^2/tau);
    sigma = sqrt((alpha + tau)/alpha/tau + sqrt(sigma));
end

s = diag(S);
lambda = s .* (s > sigma) + alpha/(alpha+tau) * s.*(s <=sigma);

if Dx > Dy
    Lambda = [diag(lambda); zeros(Dx-Dy,Dy)];
else
    Lambda = [diag(lambda) zeros(Dx,Dy-Dx)];
end

r = max(sum(lambda > 1/sqrt(tau)),1);

A = U*Lambda*V';
C = V(:,1:r) * (eye(r) - diag(1./(lambda(1:r).^2)/tau)) * V(:,1:r)';