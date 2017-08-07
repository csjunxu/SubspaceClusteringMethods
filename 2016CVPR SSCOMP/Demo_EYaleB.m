% This code tests the performance of SSC-OMP on the EYaleB face database.
% The code generates results in Table 2 of the paper
% C. You, D. Robinson, R. Vidal, Scalable Sparse Subspace Clustering by
% Orthogonal Matching Pursuit, CVPR 2016.

% In this code, we apply SSC-OMP to the face images of a randomly picked
% n \in \{2, 10, 20, 30, 38} subjects in the Extended Yale B database. Each
% subject has 64 images under different illumination conditions.

% Instructions for running the code:
% - Download code for computing clustering accuracy. Go to
% http://www.cad.zju.edu.cn/home/dengcai/Data/Clustering.html and download
% bestMap.m and Hungarian.m. Alternatively, you can use your own function
% by redefining the function evalAccuracy.m
% - Run. You can modify the parameter "nCluster" below to run for different
% number of subjects.

% Copyright Chong You @ Johns Hopkins University, 2016
% chong.you1987@gmail.com

%% Settings

%% Load data 
% Data is preprocessed and saved in the .mat file.
% EYALEB_DATA is a D by N matrix. Each column is a face image and N =
% 38 subjects * 64 images/subject = 2414. Each image is downsampled from
% 192*168 to D = 48*42 = 2016. 
% EYALEB_LABEL is a 1 by N vector. Each entry is the label for the
% corresponding column in EYALEB_DATA.

dataset = 'USPS'; % YaleB_LSR   USPS

if strcmp(dataset, 'YaleB_LSR') == 1
    load 'C:\Users\csjunxu\Desktop\SC\Datasets\YaleB_Crop.mat'   % load YaleB dataset
elseif strcmp(dataset, 'USPS') == 1
    load 'C:\Users\csjunxu\Desktop\SC\Datasets\USPS_Crop.mat'   % load USPS dataset
end

SegmentationMethod = 'SSC_OMP';
writefilepath = ['C:/Users/csjunxu/Desktop/SC/Results/' dataset '/'];

Repeat = 1; %number of repeations
DR = 1;
dim = 6;
for K = [3 4 5 6]
    for thr = [1e-8 1e-7 1e-6 1e-5 1e-4 1e-3]
        for nSet = [2 3 5 8 10];
            n = nSet;
            index = Ind{n};
            for i = 1:size(index,1)
                fea = [];
                gnd = [];
                for p = 1:n
                    fea = [fea Y{index(i, p), 1}];
                    gnd= [gnd p * ones(1, length(S{index(i, p)}))];
                end
                [D, N] = size(fea);
                fprintf( '%d: %d\n', size(index, 1), i ) ;
                missrate = zeros(size(index, 1), Repeat) ;
                for j = 1 : Repeat
                    N = length(gnd);
                    % clustering
                    tic;
                    %     fprintf('Dimension reduction...\n')
                    if DR==1
                        fea = dimReduction_PCA(fea, dim*n);
                    end
                    % normalization
                    %     fprintf('Normalization...\n')
                    fea = cnormalize_inplace(fea);
                    % generate representation
                    %     fprintf('Representation...\n')
                    R = OMP_mat_func(fea, K, thr);
                    % generate affinity
                    %     fprintf('Affinity...\n')
                    R(1:N+1:end) = 0;
                    % R = cnormalize(R, Inf);
                    A = abs(R) + abs(R)';
                    % generate label
                    %     fprintf('Generate label...\n')
                    groups = SpectralClustering(A, n, 'Eig_Solver', 'eigs');
                    time = toc;
                    % Evaluation
                    perc = evalSSR_perc( R, gnd );
                    ssr = evalSSR_error( R, gnd );
                    conn = evalConn( A, gnd);
                    accr  = evalAccuracy(gnd, groups);
                    % record
                    missrate(i, j) = 1 - accr;
                    fprintf('%.3f%% \n' , missrate(i, j)*100) ;
                end
                % output
                missrateTot{n}(i) = mean(missrate(i, :)*100);
                fprintf('Mean missrate of %d/%d is %.3f%%.\n ' , i, size(index, 1), missrateTot{n}(i)) ;
            end
            %% output
            avgmissrate(n) = mean(missrateTot{n});
            medmissrate(n) = median(missrateTot{n});
            fprintf('Total mean missrate  is %.3f%%.\n ' , avgmissrate(n)) ;
            matname = sprintf([writefilepath dataset '_' SegmentationMethod '_DR' num2str(DR) '_dim' num2str(dim) '_K' num2str(K) '_thr' num2str(thr) '.mat']);
            save(matname,'missrateTot','avgmissrate','medmissrate');
        end
    end
end
