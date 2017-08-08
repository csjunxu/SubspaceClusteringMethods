function [idx,m2] = K_flats(data, dim, idx)

% function [idx,m2] = K_flats(data, dim, idx)
%   applies the iterative K-flats algorithm to cluster data into
%   K (=length(dim)) flats of dimensions (specified by dim),
%   possibly based on an initial guess of the labels (idx).
%
% Input:
%   data: N-by-D matrix
%   dim: dimensions of the planes;
%       can be a single number if all dimensions are the same AND idx is given 
%   idx: initial labels of the points;
%       if not given, will randomly assign points to subspaces        
%
% Output:
%   idx: labels of the data points associated with the clusters
%   m2: averaged L2 error of the final model

if nargin<2
    error('Dimensions need to be specified!')
else
    K = length(dim);
end

N = size(data,1);

if nargin<3
    
    idx = ceil(K*rand(N,1));
    
%     idx = zeros(N,1);
%     seeds = randsample(N,K);
%     for k = 1:K
%         dists = sum((data-repmat(data(seeds(k),:),N,1)).^2,2);
%         [min_dist, I] = sort(dists);
%         idx(I(1:dim(k)+2)) = k*ones(dim(k)+2,1);
%     end
    
elseif K == 1
    
    K = max(idx); 

end

[ctr,dir] = computing_centers_and_bases(data,idx,dim);

m1 = 0;
m2 = L2error(data, dim, idx, ctr, dir);

if K == 1; return; end

max_loops = 1000;
tol = 1e-6;

cnt = 1;
while cnt<max_loops && abs(m1-m2)>tol

    inds_in = (idx>0);
    dist = computing_point_to_flat_distances(data(inds_in,:),ctr,dir);
    [min_dist, idx(inds_in)] = min(dist,[],2);
    
    [ctr,dir] = computing_centers_and_bases(data,idx,dim);
    
    m1 = m2;
    m2 = computing_average_L2_error(data, dim, idx, ctr, dir);

    cnt = cnt+1;

end
