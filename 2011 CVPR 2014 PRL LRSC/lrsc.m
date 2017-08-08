%--------------------------------------------------------------------------
% This is the function to call the sparse optimization program, to call the
% spectral clustering algorithm and to compute the clustering error.
% r = projection dimension, if r = 0, then no projection
% affine = use the affine constraint if true
% s = clustering ground-truth
% missrate = clustering error
% CMat = coefficient matrix obtained by SSC
%--------------------------------------------------------------------------
% Copyright @ Jun Xu, 2012
%--------------------------------------------------------------------------

function [LRSC_missrate,C] = lrsc(X,tau,r,outlier,rho,s)

if (nargin < 5)
    rho = 1;
end
if (nargin < 4)
    outlier = false;
end
if (nargin < 3)
    r = 0;
end
if (nargin < 2)
    tau = 100/norm(D)^2;
    alpha = 0.5*tau;
end

n = max(s);
Xp = DataProjection(X,r);

if (~outlier)
    %Cp = lrsc_noiseless(Xp,tau);
    [A,Cp] = lrsc_noisy(Xp,tau);
else
    [Cp,A,E] = lrsc_outliersrelax(Xp,tau,alpha,gamma);
end

C = BuildAdjacency(thrC(Cp,rho));
groups = SpectralClustering(C,n);
%compute misclassification rate
LRSC_missrate = Misclassification(groups,s);