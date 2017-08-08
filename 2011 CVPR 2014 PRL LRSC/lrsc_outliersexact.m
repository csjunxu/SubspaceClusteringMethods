%--------------------------------------------------------------------------
% [C,A,E] = lrsc_outliersexact(D,tau,gamma,Nit,whichNorm)
% Low Rank Subspace Clustering algorithm for data lying in a union of 
% subspaces and contaminated with outliers
%
% min |C|_* + tau/2*|A-AC|_F^2 + gamma*|E|_1 s.t. C = C', D = A + E
%
% C = affinity matrix
% A = clean data matrix whose columns are points in a union of subspaces
% E = matrix of sparse outliers
% D = noisy data matrix whose columns are points in a union of subspaces
% tau = scalar parameter 
% gamma = scalar parameter
%--------------------------------------------------------------------------
% Copyright @ Paolo Favaro, December 2012
%--------------------------------------------------------------------------

function [C,A,E] = lrsc_outliersexact(D,tau,gamma,Nit,whichNorm)

if nargin<5
    whichNorm = 'L1';
end
if nargin<4
    Nit = 150;
end
if Nit<1
    Nit = 150;
end


% initialization
E = zeros(size(D));
Y = zeros(size(D));
mu = 1e2; % used for the L2-norm augmented Lagrangian term on D=A+E
rho = 1.1;
for it = 1:Nit

    %%% COMPUTE A
    A = lrsc_noisy(D-E+Y/mu,tau,mu);
    
    %%% COMPUTE outliers E
    [E,~] = mregularize(D-A+Y/mu,whichNorm,gamma/mu);
    
    %%% COMPUTE Lagrange multipliers Y
    Yprev = Y;
    Y = Y+mu*(D-A-E);
    
    mu = mu*rho;
    
    %%% stop if no significant change occurred
    if mean(abs(Y(:)-Yprev(:)))<1e-7
        figure(3)
        hold all
        plot(it,rand(1),'o');
        break;
    end
    
end
[A,C] = lrsc_noisy(D-E+Y/mu,tau,mu);