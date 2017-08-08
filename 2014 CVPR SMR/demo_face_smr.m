clear;
warning off all

load 'C:\Users\csjunxu\Desktop\SC\Datasets\YaleB_Crop.mat'              % load YaleB dataset

dataset = 'YaleB_LSR';
writefilepath = ['C:/Users/csjunxu/Desktop/SC/Results/' dataset '/'];

Repeat = 1; %number of repeations
DR = 1; % perform dimension reduction or not
if DR == 0
    dim = size(Y{1, 1}, 1);
elseif DR == 1
    dim = 6;
else
    DR = 1;
    dim = 6;
end

SegmentationMethod = 'SMR' ;

% para.knn = 4;
% para.gamma = 5;
% para.elpson = 0.01;
% para.aff_type = 'J1';
% para.alpha = 2^15;

%% Subspace segmentation
Par.aff_type = 'J1';
for gamma = [5]
    Par.gamma = gamma;
    for knn = [4]
        Par.knn = knn;
        for elpson = [2]
            Par.elpson = 10^(-elpson);
            for alpha = [15]
                Par.alpha = 2^(alpha);
                for nSet = [2 3 5 8 10]
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
                        redDim = size(fea, 1);
                        if DR == 1
                            %% PCA Projection
                            [ eigvector , eigvalue ] = PCA( fea ) ;
                            maxDim = length(eigvalue);
                            fea = eigvector' * fea ;
                            redDim = min(n*dim, size(fea, 1)) ;
                        end
                        %% normalize
                        for c = 1 : size(fea,2)
                            fea(:,c) = fea(:,c) /norm(fea(:,c)) ;
                        end
                        %% Subspace Clustering
                        missrate = zeros(size(index, 1), Repeat) ;
                        fprintf( 'dimension = %d \n', redDim ) ;
                        Yfea = fea(1:redDim, :) ;
                        for j = 1 : Repeat
                            C = smr( Yfea , Par ) ;
                            for k = 1 : size(C,2)
                                C(:, k) = C(:, k) / max(abs(C(:, k))) ;
                            end
                            Z = ( abs(C) + abs(C') ) / 2 ;
                            idx = clu_ncut(Z,n) ;
                            missrate(i, j) = 1 - compacc(idx,gnd);
                            fprintf('%.3f%% \n' , missrate(i, j)*100) ;
                        end
                        missrateTot{n}(i) = mean(missrate(i, :)*100);
                        fprintf('Mean error of %d/%d is %.3f%%.\n ' , i, size(index, 1), missrateTot{n}(i)) ;
                    end
                    %% output
                    avgmissrate(n) = mean(missrateTot{n});
                    medmissrate(n) = median(missrateTot{n});
                    fprintf('Total mean error  is %.3f%%.\n ' , avgmissrate(n)) ;
                    allavgmissrate = mean(avgmissrate(avgmissrate~=0));
                    matname = sprintf([writefilepath dataset '_' SegmentationMethod '_knn' num2str(Par.knn) '_gamma' num2str(Par.gamma) '_elpson' num2str(Par.elpson) '_aff' Par.aff_type '_alpha' num2str(Par.alpha) '.mat']);
                    save(matname,'avgmissrate','medmissrate', 'allavgmissrate');
                end
            end
        end
    end
end


