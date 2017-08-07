function [perc, x] = evalSSR_perc(C, s)
% C \in R^N-by-N: representation matrix by self-expressiveness based method
% s \in {1, 2, ... n}^N: group labels
%     
% perc: percentage of points that is subspace sparse.
% x: for plotting curve: [0; x; 1] \times [0:1/N:1, 1];

N = length(s);
x = zeros(N, 1);
for ii = 1:N
    x(ii) = max( abs(C(s ~= s(ii), ii)) );
end
x = x ./ max(abs(C), [], 1)';
perc = sum(x < 1e-5) / N;
