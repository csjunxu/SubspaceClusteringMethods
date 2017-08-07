clear ;

%% reduced dimension
ProjRank = 12 ;
datadir = 'C:/Users/csjunxu/Desktop/SC/Datasets/Hopkins155/';
seqs = dir(datadir);
% Get rid of the two directories: "." and ".."
seq3 = seqs(3:end);
% Save the data loaded in struct "data "
data = struct('ProjX', {}, 'name',{}, 'ids',{});

dataset = 'Hopkins155';

resultdir = 'C:/Users/csjunxu/Desktop/SC/Results/';

for i=1:length(seq3)
    fname = seq3(i).name;
    fdir = [datadir '/' fname];
    if isdir(fdir)
        datai = load([fdir '/' fname '_truth.mat']);
        id = length(data)+1;
        % the true group numbers
        data(id).ids = datai.s;
        % file name
        data(id).name = lower(fname);
        % X is the motion sequence
        X = reshape(permute(datai.x(1:2,:,:),[1 3 2]), 2*datai.frames, datai.points);
        
        % PCA projection
        [ eigvector , eigvalue ] = PCA( X ) ;
        ProjX = eigvector(:,1:ProjRank)' * X ;
        data(id).ProjX = [ProjX ; ones(1,size(ProjX,2)) ] ;
    end
end
clear seq3;


%% Subspace segmentation methods

SegmentationMethod = 'SSC_OMP' ;

for K = [4 5 6 7 8 9 10]
    for thr = [1e-8 1e-7 1e-6 1e-5 1e-4 1e-3]
        maxNumGroup = 5;
        for i = 1:maxNumGroup
            num(i) = 0;
        end
        %%
        errs = zeros(length(data),1);
        for i = 1 : length(data)
            ProjX = data(i).ProjX ;
            [D, N] = size(ProjX);
            gnd = data(i).ids' ;
            K = length( unique( gnd ) ) ;
            n = max(gnd);
            C = OMP_mat_func(ProjX, K, thr);
            % generate affinity
            %     fprintf('Affinity...\n')
            C(1:N+1:end) = 0;
            nCluster = length( unique( gnd ) ) ;
            Z = ( abs(C) + abs(C') ) / 2 ;
            idx = clu_ncut(Z,nCluster) ;
            accuracy = compacc(idx,gnd) ;
            missrate = 1-accuracy;
            num(n) = num(n) + 1;
            missrateTot{n}(num(n)) = missrate;
            fprintf('seq %d\t %f\n', i , missrate ) ;
        end
        fprintf('\n') ;
        
        L = [2 3];
        allmissrate = [];
        for i = 1:length(L)
            j = L(i);
            avgmissrate(j) = mean(missrateTot{j});
            medmissrate(j) = median(missrateTot{j});
            allmissrate = [allmissrate missrateTot{j}];
        end
        avgallmissrate = sum(allmissrate)/length(allmissrate);
        medallmissrate = median(allmissrate);
        matname = sprintf([resultdir dataset '_' SegmentationMethod '_K' num2str(K) '_thr' num2str(thr) '.mat']);
        save(matname,'avgallmissrate','medallmissrate','missrateTot','avgmissrate','medmissrate');
    end
end




