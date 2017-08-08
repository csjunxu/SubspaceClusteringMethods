%--------------------------------------------------------------------------
% [A,C] = polythreshAP(Delta,tau,alpha)
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
% Copyright @ Paolo Favaro, January 2013
%--------------------------------------------------------------------------
function [A,C] = polythreshAP(Delta,tau,alpha)
% compute polynomial thresholding operator numerically
% notice that the singular values of A are always smaller
% than the singular values of D

[U,S,V] = svd(Delta,'econ');
sigma = diag(S);

if max(sigma)>0
    lambdas = 0:max(sigma)/300:max(sigma);
    [lam,sig] = meshgrid(lambdas,sigma);
    cost = alpha/2*(sig-lam).^2;
    I1 = lam>1/sqrt(tau);
    if tau>1e-6
        cost(I1) = cost(I1)+(1-1/2/tau*lam(I1).^-2);
    end
    I2 = lam<=1/sqrt(tau);
    if tau<1e28
        cost(I2) = cost(I2)+tau/2*lam(I2).^2;
    end
    [~,ind] = min(cost,[],2);
    lambda = lambdas(ind);
else
    lambda = sigma;
end
A = U*diag(lambda)*V';

if nargout>1
    n = sum(lambda>1/sqrt(tau));
    if n>0
        C = V(:,1:n)*(eye(n)-1/tau*diag(1./lambda(1:n).^2))*V(:,1:n)';
    else
        C = zeros(size(Ak,2));
    end
end