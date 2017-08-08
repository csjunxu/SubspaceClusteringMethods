% @input:   X - d x n data matrix
%           k - number of blocks/subspaces
%           lambda - trade-off parameter
%           sdim - total dimension of the underlying subspace
% @output:  Z - representation coefficient matrix
%           E - noise matrix
% author: jshfeng@gmail.com

function [Z,E] = ssgd(X,k,lambda,sdim,Z0, SegMethod, ProjMethod, T, beta)
[d,n] = size(X);

%% initialization
Z = Z0;
E = zeros(d,n);
I = eye(n);
XXt = X'*X;

if strcmp(ProjMethod,'lsr')
%     alpha = 4.6e-3;
    alpha = 0.4;
    S = (XXt+alpha*I)\XXt;
    W = abs(S)+abs(S');
    idx = clu_ncut(W,k);
else
    idx = [];
    S = [];
end

%% parameters for learning rate
beta = 1;
T = 1000;
p = sdim*2;
r = sdim;
delta = 1.5*n;
sigma = norm(XXt);
eta = beta*sqrt(p)*delta;
gamma = delta;

switch SegMethod
    case 'lrr'
        eta = eta/(sqrt(n)*(sigma*gamma*lambda+sqrt(r))*sqrt(T));
%             eta = 0.005;
    case 'ssc'
        eta = eta/(sqrt(n)*(sigma*gamma*lambda+n))*sqrt(T);
    case 'lsr'
        eta = eta/(sqrt(n)*(sigma*gamma*lambda+gamma))*sqrt(T);
end

%% optimization
isRandomized = 0;

for iter = 1:T
    
    %% sub-gradient calculation
    switch SegMethod
        case 'lrr'
            [U,junk,V] = svd(Z,'econ');
            grad = lambda*XXt*(Z-I) + U*V'; % subgradient
        case 'ssc'
            Zs = sign(Z);
%             Zs(Z==0) = 1;
            grad = lambda*XXt*(Z-I) + Zs;
            grad = grad - diag(diag(grad));
        case 'lsr'
            grad = lambda*XXt*(Z-I) + 2*Z;
    end
    
    if isRandomized
        ind = randperm(n,p);
        Y = I(:,ind);
        Y = sqrt(n)/sqrt(p)*Y;
        grad = grad*Y*Y'; % randomized subgradient
    end

    %% gradient descent
    Z = Z - eta*grad;
%     Z = projKappa(Z,k,idx,S,ProjMethod,beta); % projection
end
end

%% projection function
function Z = projKappa(Z0,k,idx,S,method,beta) 

Z = zeros(size(Z0));

if strcmp(method,'ncut')
    W = 0.5*(abs(Z0) + abs(Z0'));
    NcutDiscrete = ncutW(W,k);
end

if strcmp(method,'lsr')
    temp = Z0 + beta*S;
    W = 0.5*(abs(temp) + abs(temp'));
    NcutDiscrete = ncutW(W,k);
end

for i = 1:k
    if strcmp(method,'ncut')
        index = find(NcutDiscrete(:,i));
    else
%         index = find(idx==i);
        index = find(NcutDiscrete(:,i));
    end
    Z(index,index) = Z0(index,index);
end
end