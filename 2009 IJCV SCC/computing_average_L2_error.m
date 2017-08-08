function mse = computing_average_L2_error(data, dim, idx, ctr, dir)

N = size(data,1);

K = max(idx);

if length(dim) < K
    dim = dim(end)*ones(K,1);
end

if nargin<4
    [ctr,dir] = computing_centers_and_bases(data,idx,dim);
end

mse = 0;
for k = 1:K
    cls_k = data((idx==k),:);
    n_k = size(cls_k,1);
    if n_k > dim(k)
        cls_k = cls_k - repmat(ctr{k,1},n_k,1);
        mse = mse + sum(sum(cls_k.^2,2) - sum((cls_k*dir{k,1}').^2,2));
    end
end

mse = sqrt(mse/N);