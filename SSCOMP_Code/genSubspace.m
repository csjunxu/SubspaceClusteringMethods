function [X, s] = genSubspace(D, n, Ni, di, sigma, corruption)
% This code generates a matrix X of D by sum(Ni) that contains n subspaces
% of dimension di, with noise level sigma.
% Input:
%     D: dimension of ambience space
%     n: number of subspace
%     Ni: #points in each subspace
%     di: dimension of each subspace
%     sigma: noise deviation
% Output:
%     X: data (D by sum(Ni))
%     s: label of X (1 by sum(Ni))
% example: 
%     [X, s] = genSubspace(100, 6, 100, 5, 0);

% Copyright Chong You @ Johns Hopkins University, 2016
% chong.you1987@gmail.com

if length(Ni) == 1
    Ni = repmat(Ni, 1, n);
end
if length(di) == 1
    di = repmat(di, 1, n);
end
if ~exist('sigma', 'var')
    sigma = 0;
end
if ~exist('corruption', 'var')
    corruption = 0;
end

X = zeros(D, sum(Ni)); s = zeros(1, sum(Ni));
idx = 0;
for in = 1:n
%     Xtmp = randn(D, di(in)) * randn(di(in), Ni(in));

    Xtmp = randn(D, D);
    [Utmp, ~, ~] = svds(Xtmp, di(in)); % generate a random subspace
    Vtmp = randn(di(in), Ni(in));
    Xtmp = Utmp * Vtmp; % generate random points in subspace
    Xtmp = bsxfun(@rdivide, Xtmp, sqrt(sum(Xtmp .^2, 1)) ); % normalize
    
    X(:, [idx+1 : idx+Ni(in)]) = Xtmp;
    s([idx+1 : idx+Ni(in)]) = in;
    idx = idx+Ni(in);
end

noise_term = sigma * randn(D, sum(Ni)) / sqrt(D);
X = X + noise_term;
corruption_mask = randperm( D*sum(Ni), round( corruption*D*sum(Ni) ) );
X(corruption_mask) = 0;

end