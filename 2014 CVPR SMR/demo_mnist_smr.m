clear;

addpath('MNISThelpcode');
addpath('C:\Users\csjunxu\Documents\GitHub\SubspaceCluteringCode\SSCOMP_Code\scatnet-0.2');
%% Settings
for nSample = [200] % number of images for each digit
    
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
            MNIST_DATA = loadMNISTImages('train-images.idx3-ubyte');
            MNIST_LABEL = loadMNISTLabels('train-labels.idx1-ubyte');
            MNIST_SC_DATA = SCofDigits(MNIST_DATA);
            save C:\Users\csjunxu\Desktop\SC\Datasets\MNIST_SC.mat MNIST_SC_DATA MNIST_LABEL;
        end
        MNIST_DATA = MNIST_SC_DATA;
    end
    
    dataset = 'MNIST';
    writefilepath = ['C:/Users/csjunxu/Desktop/SC/Results/' dataset '/'];
    
    nExperiment = 20; % number of repeations
    DR = 1; % perform dimension reduction or not
    if DR == 0
        dim = size(Y{1, 1}, 1);
    elseif DR == 1
        dim = 50;
    else
        DR = 1;
        dim = 50;
    end
    
    %% Subspace segmentation methods
    SegmentationMethod = 'SMR' ;
    
    %% Subspace segmentation
    Par.aff_type = 'J1';
    Par.gamma = 1;
    for knn = [4]
        Par.knn = knn;
        for elpson = [-3]
            Par.elpson = 10.^elpson;
            for alpha = [-7 -8 -9]
                Par.alpha = 2.^alpha;
                missrate = zeros(nExperiment, 1) ;
                for i = 1:nExperiment
                    nCluster = 10;
                    digit_set = 0:9; % set of digits to test on, e.g. [2, 0]. Pick randomly if empty.
                    % prepare data
                    if isempty(digit_set)
                        rng(i); Digits = randperm(10, nCluster) - 1;
                    else
                        Digits = digit_set;
                    end
                    if length(nSample) == 1
                        nSample = ones(1, nCluster) * nSample;
                    end
                    mask = zeros(1, sum(nSample));
                    gnd = zeros(1, sum(nSample));
                    nSample_cum = [0, cumsum(nSample)];
                    for iK = 1:nCluster % randomly take data for each digit.
                        allpos = find( MNIST_LABEL == Digits(iK) );
                        rng( (i-1) * nCluster + iK );
                        selpos = allpos( randperm(length(allpos), nSample(iK)) );
                        
                        mask( nSample_cum(iK) + 1 : nSample_cum(iK+1) ) = selpos;
                        gnd( nSample_cum(iK) + 1 : nSample_cum(iK+1) ) = iK * ones(1, nSample(iK));
                    end
                    fea = MNIST_DATA(:, mask);
                    N = length(gnd);
                    
                    %% PCA Projection
                    redDim = size(fea, 1);
                    if DR == 1
                        [ eigvector , eigvalue ] = PCA( fea ) ;
                        maxDim = length(eigvalue) ;
                        fea = eigvector' * fea ;
                        redDim = min(nCluster*dim, size(fea, 1)) ;
                    end
                    %                     fprintf( 'dimension = %d \n', redDim ) ;
                    %% normalize
                    for c = 1 : size(fea,2)
                        fea(:,c) = fea(:,c) /norm(fea(:,c)) ;
                    end
                    %% Subspace Clustering
                    Yfea = fea(1:redDim, :) ;
                    C = smr(Yfea, Par);
                    %% generate affinity
                    for k = 1 : size(C, 2)
                        C(:, k) = C(:, k) / max(abs(C(:, k))) ;
                    end
                    Z = ( abs(C) + abs(C') ) / 2 ; % abs is useless in our model
                    %% generate label
                    idx = clu_ncut(Z, nCluster) ;
                    %% Evaluation
                    missrate(i) = 1 - compacc(idx, gnd);
                    fprintf('%d: %.3f%% \n' , i, missrate(i)*100) ;
                end
                %% output
                avgmissrate = mean(missrate*100);
                medmissrate = median(missrate*100);
                fprintf('Total mean missrate  is %.3f%%.\n' , avgmissrate) ;
                matname = sprintf([writefilepath dataset '_' num2str(nSample(1)) '_' num2str(nExperiment) '_' SegmentationMethod '_knn' num2str(Par.knn) '_elpson10^' num2str(elpson) '_alpha2^' num2str(alpha) '_aff' Par.aff_type '_gamma' num2str(Par.gamma) '.mat']);
                save(matname,'missrate','avgmissrate','medmissrate');
            end
        end
    end
end
