% This code tests the performance of SSC-OMP on the MNIST digit database.
% The code generates results in Table 1 of the paper
% C. You, D. Robinson, R. Vidal, Scalable Sparse Subspace Clustering by 
% Orthogonal Matching Pursuit, CVPR 2016.

% In this code, we apply SSC-OMP to the handwritten digit images. We 
% randomly pick N_i \in \{50, 100, 200, 400, 600\} digits from each of the
% 10 digits (i.e. 0-9), with feature extracted from a scattering network
% and projected to dimension 500.

% Instructions for running the code:
% - Download code for computing clustering accuracy. Go to
% http://www.cad.zju.edu.cn/home/dengcai/Data/Clustering.html and download
% bestMap.m and Hungarian.m. Alternatively, you can use your own function
% by redefining the function evalAccuracy.m
% - Download data. We use all 60,000 training images from MNIST and apply 
% scattering transform. Specifically, 
%   = Download MNIST traning image/label files from 
%     http://yann.lecun.com/exdb/mnist/ and put the data in folder MNIST/
%   = Download files for reading data (loadMNISTImages.m and 
%     loadMNISTLabels.m) from
%     http://ufldl.stanford.edu/wiki/index.php/MATLAB_Modules
%   = Download scattering transform package ScatNet (v0.2) from
%     http://www.di.ens.fr/data/software/
%     and install.
% - Run. You can modify the parameter "nSample" below to run for different
% number of samples per digit. 

% Copyright Chong You @ Johns Hopkins University, 2016
% chong.you1987@gmail.com

%% Settings
% setup
nSample = 400; % number of images for each digit

% dimension reduction
reduceDimension = @(data) dimReduction_PCA(data, 500);
% normalization
normalizeColumn = @(data) data;
% representation     
buildRepresentation = @(data) OMP_mat_func(data, 10, 1e-6); % second parameter is sparsity
% spectral clustering       
genLabel = @(affinity, nCluster) SpectralClustering(affinity, nCluster, 'Eig_Solver', 'eigs');

%% Load data
addpath('C:\Users\csjunxu\Desktop\SC\Datasets\MNIST\')
if ~exist('MNIST_DATA', 'var')
    try
        % MNIST_SC_DATA is a D by N matrix. Each column contains a feature 
        % vector of a digit image and N = 60,000.
        % MNIST_LABEL is a 1 by N vector. Each entry is the label for the
        % corresponding column in MNIST_SC_DATA.
        load MNIST_SC.mat MNIST_SC_DATA MNIST_LABEL;
    catch
        MNIST_DATA = loadMNISTImages('train-images-idx3-ubyte');
        MNIST_LABEL = loadMNISTLabels('train-labels-idx1-ubyte'); 
        MNIST_SC_DATA = SCofDigits(MNIST_DATA);
        save MNIST_SC.mat MNIST_SC_DATA MNIST_LABEL;
    end
    MNIST_DATA = MNIST_SC_DATA;
end

%% Clustering
nExperiment = 20;
results = zeros(nExperiment, 6); %results
for iExperiment = 1:nExperiment
    nCluster = 10;
    digit_set = 0:9; % set of digits to test on, e.g. [2, 0]. Pick randomly if empty.
    % prepare data
    if isempty(digit_set)
        rng(iExperiment); Digits = randperm(10, nCluster) - 1;
    else
        Digits = digit_set;
    end
    if length(nSample) == 1
        nSample = ones(1, nCluster) * nSample;
    end
    mask = zeros(1, sum(nSample));
    s = zeros(1, sum(nSample));
    nSample_cum = [0, cumsum(nSample)];
    for iK = 1:nCluster % randomly take data for each digit.
        allpos = find( MNIST_LABEL == Digits(iK) );
        rng( (iExperiment-1) * nCluster + iK );
        selpos = allpos( randperm(length(allpos), nSample(iK)) );

        mask( nSample_cum(iK) + 1 : nSample_cum(iK+1) ) = selpos;
        s( nSample_cum(iK) + 1 : nSample_cum(iK+1) ) = iK * ones(1, nSample(iK));
    end
    X = MNIST_DATA(:, mask);
    N = length(s);
    
    % Clustering
    tic;
    
    fprintf('Dimension reduction...\n')
    X = reduceDimension(X);
    % normalization
    fprintf('Normalization...\n')
    X = normalizeColumn(X);
    % generate representation
    fprintf('Representation...\n')
    R = buildRepresentation(X);
    % generate affinity
    fprintf('Affinity...\n')
    R(1:N+1:end) = 0;
    R = cnormalize(R, Inf);
    A = abs(R) + abs(R)';
    % generate label
    fprintf('Generate label...\n')
    groups = genLabel(A, nCluster);                             

    time = toc;
    
    % Evaluation
    perc = evalSSR_perc( R, s );
    ssr = evalSSR_error( R, s );   
    conn = evalConn( A, s);
    accr  = evalAccuracy(s, groups);

    % output
    dataformat = '%d-th experiment: perc = %f, ssr = %f, conn = %f, accr = %f, time = %f\n';
    dataValue = [iExperiment, perc, ssr, conn, accr, time];
    fprintf(dataformat, dataValue);
    % record
    results(iExperiment, :) = dataValue;
end
% output
dataValue = mean(results, 1);
fprintf('\nAverage: perc = %f, ssr = %f, conn = %f, accr = %f, time = %f\n', dataValue(2:end));

