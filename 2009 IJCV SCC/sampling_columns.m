function inds = sampling_columns(sampleLabels,K,opts)

%K = max(sampleLabels); % number of clusters
%inliers = find(sampleLabels>0); % remove outliers

c_ave = opts.c/K; 

inds = zeros(opts.n-1,opts.c);
for k = 1:K

    % indices of points in cluster k and its size
    inds_k = find(sampleLabels == k);
    n_k = length(inds_k);

    if n_k < opts.n-1 % not enough points in cluster k
        inds_k = find(sampleLabels>0);
    end
    
    % sample c_ave columns from cluster k
    for col = ((k-1)*c_ave+1):(k*c_ave) 
        inds(:,col) = randsample(inds_k,opts.n-1);
    end

end

% for col = (K*c_ave+1):opts.c
%     inds(:,col) = randsample(inliers,opts.n-1);
% end