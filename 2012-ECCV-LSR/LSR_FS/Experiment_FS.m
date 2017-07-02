clear ;
addpath('LSR');
load 'C:\Users\csjunxu\Desktop\SC\Datasets\YaleB_Crop.mat'              % load YaleB dataset
% load 'C:\Users\csjunxu\Desktop\SC\Datasets\USPS_Crop.mat'   % load USPS dataset

writefilepath = ['C:/Users/csjunxu/Desktop/SC/Results/' dataset '/'];

%% Subspace segmentation methods
% SegmentationMethod = 'LSR1' ;     % LSR1 by (16) in our paper
SegmentationMethod = 'LSR2' ;     % LSR2 by (18) in our paper

Repeat = 1; %number of repeations
DR = 1; % perform dimension reduction or not
dim = 6;

%% Subspace segmentation
%% Parameter
% switch n
%     case 5
%         para = [0.4] * ones(1,Repeat) ;
%     case 10
%         para = [0.004 ] * ones(1,Repeat) ;
% end

for lambda = [.05 .06 .04]
    para = lambda * ones(1,Repeat) ;
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
            redDim = size(fea, 1);
            if DR == 1
                %% PCA Projection
                [ eigvector , eigvalue ] = PCA( fea ) ;
                maxDim = length(eigvalue);
                fea = eigvector' * fea ;
                redDim = n * dim ;
            end
            %% normalize
            for c = 1 : size(fea,2)
                fea(:,c) = fea(:,c) /norm(fea(:,c)) ;
            end
            
            %% Subspace segmentation
            missrate = zeros(size(index, 1), Repeat) ;
            
            fprintf( 'dimension = %d \n', redDim ) ;
            Yfea = fea(1:redDim,:) ;
            for j = 1 : Repeat
                p = para( j ) ;
                Accuracy(i,j) = SubspaceSegmentation( SegmentationMethod , Yfea , gnd , p ) ;
                missrate(i, j) = 1 - Accuracy(i,j);
                fprintf('%.3f%% \n' , missrate(i, j)*100) ;
            end
            
            missrateTot{n}(i) = mean(missrate(i, :)*100);
            fprintf('Mean Accuracy of %d/%d is %.3f%%.\n ' , i, size(index, 1), missrateTot{n}(i)) ;
        end
        %% output
        avgmissrate(n) = mean(missrateTot{n});
        medmissrate(n) = median(missrateTot{n});
        fprintf('Total mean missrate  is %.3f%%.\n ' , avgmissrate(n)) ;
        matname = sprintf([writefilepath 'YaleB_Crop_' SegmentationMethod '_DR' num2str(DR) '_dim' num2str(dim) '_lambda' num2str(lambda) '.mat']);
        save(matname,'missrateTot','avgmissrate','medmissrate');
    end
end
