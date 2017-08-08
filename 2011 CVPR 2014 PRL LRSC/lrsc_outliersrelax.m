%--------------------------------------------------------------------------
% [C,A,E] = lrsc_outliersrelax(D,tau,alpha,gamma,Nit,whichNorm)
% Low Rank Subspace Clustering algorithm for data lying in a union of 
% subspaces and contaminated with outliers
%
% min |C|_* + tau/2*|A-AC|_F^2 + alpha/2*|D-A-E|_F^2 + gamma*|E|_1 s.t. C = C'
%
% C = affinity matrix
% A = clean data matrix whose columns are points in a union of subspaces
% E = matrix of sparse outliers
% D = noisy data matrix whose columns are points in a union of subspaces
% tau = scalar parameter 
% alpha = scalar parameter 
% gamma = scalar parameter
%--------------------------------------------------------------------------
% Copyright @ Paolo Favaro, January 2013
%--------------------------------------------------------------------------

function [C,A,E] = lrsc_outliersrelax(D,tau,alpha,gamma,Nit,whichNorm)

if nargin<5
    Nit = 150;
end
if nargin < 4
    gamma = 1; % not defined
end
if nargin < 3
    alpha = 0.5*tau;
end
if nargin < 2
    tau = 100/norm(D)^2;
    alpha = 0.5*tau;
end

if Nit<1
    Nit = 150;
end

if nargin<6
    whichNorm = 'L1';
end

% initialization
E = zeros(size(D));
Aprev = D;
for it = 1:Nit

    %%% COMPUTE A
    A = lrsc_noisy(D-E,tau,alpha);
    
    %%% COMPUTE outliers E
    [E,~] = mregularize(D-A,whichNorm,gamma/alpha);
    
    %%% stop if no significant change occurred
    if mean(abs(A(:)-Aprev(:)))<1e-7
%         figure(3)
%         hold all
%         plot(it,rand(1),'o');
        break;
    end
    Aprev = A;
    
end
[A,C] = lrsc_noisy(D-E,tau,alpha);
