function [ ids ] = spectral_clustering(W, k)

D = diag(1./sqrt(sum(W, 2)));
W = D * W * D;
[U, S, V] = svd(W);
V = U(:, 1 : k);
V = normr(V);

ids = kmeans(V, k, 'emptyaction', 'singleton', 'replicates', 10, 'display', 'off');
% ids = kmeans(V, k, 'start','sample','maxiter', 1000,'replicates',100,'EmptyAction','singleton');
end
