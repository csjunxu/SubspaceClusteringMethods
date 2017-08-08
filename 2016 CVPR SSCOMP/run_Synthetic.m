% This code tests the performance of SSC-OMP on the synthetic database.
% The code generates results in Figure 1 of the paper
% C. You, D. Robinson, R. Vidal, Scalable Sparse Subspace Clustering by 
% Orthogonal Matching Pursuit, CVPR 2016.

% The data are drawn from 5 subspaces of dimension 6 in ambient dimension
% 9. Each subspace contains the same number of points and the overall 
% number of points is varied from 150 to 10^5.

% Instructions for running the code:
% - Download code for computing clustering accuracy. Go to
% http://www.cad.zju.edu.cn/home/dengcai/Data/Clustering.html and download
% bestMap.m and Hungarian.m. Alternatively, you can use your own function
% by redefining the function evalAccuracy.m
% - Run. 

% Copyright Chong You @ Johns Hopkins University, 2016
% chong.you1987@gmail.com

%% Setting

% data
D = 9; % ambient dimension
n = 5; % number of subspaces
di = 6; % dimension of subspaces
% Ni = 50; % number of points per subspace
sigma = 0.00; % noise level

nExperiment = 20;

buildRepresentation = @(data) OMP_mat_func(data, 6, 1e-6); % second parameter is sparsity
genLabel = @(affinity, nCluster) SpectralClustering(affinity, nCluster, 'Eig_Solver', 'eigs');

%% Test begin
Ni_list = logspace(log10(30), log10(20000), 12);
for ii = 1:length(Ni_list)
    Ni = ceil(Ni_list(ii));

    results = zeros(nExperiment, 6);
    for iExperiment = 1:nExperiment
        % Generate subspaces
        rng(iExperiment);
        if length(Ni) == 1
            Ni = repmat(Ni, 1, n);
        end
        [X, s] = genSubspace(D, n, Ni, di, sigma);
        N = sum(Ni);

        % Clustering
        tic; 
        R = buildRepresentation(X);
        time1 = toc;
        R(1:N+1:end) = 0;
        A = abs(R) + abs(R)';
        groups = genLabel(A, n);             
        time2 = toc;
        
        % Evaluation
        [perc, vec] = evalSSR_perc( R, s );
        ssr = evalSSR_error( R, s );   
        conn = evalConn( A, s);
        accr  = evalAccuracy(s, groups);
        
        % output
%         dataformat = '%d-th experiment: perc = %f, ssr = %f, conn = %f, accr = %f, time1 = %f, time2 = %f\n';
        dataValue = [iExperiment, perc, ssr, conn, accr, time1, time2];
%         fprintf(dataformat, dataValue);
%         fprintf('\n');

        % record
        results(iExperiment, :) = dataValue(2:end);
    end

    % output to file
    dataformat = 'Ni = %d: perc = %f, ssr = %f, conn = %f, accr = %f, time1 = %f, time2 = %f\n';
    dataValue = [Ni(1), mean(results, 1)];
    fprintf(dataformat, dataValue);
    fprintf('\n');

end
