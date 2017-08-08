clc;
addpath('C:\Users\csjunxu\Desktop\SC\2009 CVPR 2013 PAMI SSC');
addpath('C:\Users\csjunxu\Desktop\SC\2010 ICML 2013 PAMI LRR\code\');
addpath('C:\Users\csjunxu\Documents\GitHub\Non-negativeSubspaceClustering');

%% data pre-processing
datadir = 'C:/Users/csjunxu/Desktop/SC/Datasets/Hopkins155/';
if ~exist([datadir,'Hopkins155.mat'],'file')
    seqs = dir(datadir);
    seq3 = seqs(3:end);
    data = struct('X',{},'name',{},'ids',{});
    for i=1:length(seq3)
        fname = seq3(i).name;
        fdir = [datadir '/' fname];
        if isdir(fdir)
            datai = load([fdir '/' fname '_truth.mat']);
            id = length(data)+1;
            data(id).ids = datai.s;
            data(id).name = lower(fname);
            X = reshape(permute(datai.x(1:2,:,:),[1 3 2]),...
                2*datai.frames,datai.points);
            data(id).X = [X;ones(1,size(X,2))];
        end
    end
    clear seq3;
    save([datadir,'Hopkins155.mat'],'data','-v7.3');
else
    load([datadir,'Hopkins155.mat'],'data');
end

%% parameter setting

SegMethod = 'lrr'; % the method for subspace segmentation
redo = 0;
plusLC = 0;
postProc = 0;

% parameters for LC-BD
ProjMethod = 'ncut'; % the method for the projection step
lambdaLC = 10;

% parameters for LRR
lambda = 2.4;

% parameters for SSC
alpha = 800;
r = 0;
affine = true;
outlier = false;
rho = 0.7;

% parameters  for LSR
lambda_lsr = 4e-3;
 
% path setting
dataset = 'Hopkins155';
resultPath = ['C:/Users/csjunxu/Desktop/SC/Results/' dataset '/'];


%% segmentation 
errs = zeros(length(data),1);
fprintf('start segmentation by %s... \n', SegMethod);
for i = 1:length(data)
    X = data(i).X;
    gnd = data(i).ids;
    K = max(gnd);
    if abs(K-2)>0.1 && abs(K-3)>0.1
        id = i; % the discarded sequqnce
    end
    
    switch SegMethod
        case 'lrr'            
            fname = sprintf('%s_recovered_matrix_%d.mat',SegMethod,i);
            if ~exist([resultPath,SegMethod,'MatrixZ/', fname],'file') || redo == 1
                Z = solve_lrr(X,lambda);
                save([resultPath,SegMethod,'MatrixZ/', fname], 'Z','-v7.3');
            else
                load([resultPath,SegMethod,'MatrixZ/', fname], 'Z');
            end
        case 'ssc'
            fname = sprintf('%s_recovered_matrix_%d.mat',SegMethod,i);
            if ~exist([resultPath,SegMethod,'MatrixZ/', fname],'file') || redo == 1
                Z = SSC(X,r,affine,alpha,outlier,rho,gnd);
                save([resultPath,SegMethod,'MatrixZ/',fname],'Z','-v7.3');
            else
                load([resultPath,SegMethod,'MatrixZ/',fname],'Z');
            end
        case 'lsr'
            fname = sprintf('%s_recovered_matrix_%d.mat',SegMethod,i);
            if ~exist([resultPath,SegMethod,'MatrixZ/', fname],'file') || redo == 1
                XtX = X'*X;
                Z = (XtX+lambda_lsr*eye(size(X,2)))\(XtX);
                save([resultPath,SegMethod,'MatrixZ/',fname],'Z','-v7.3');
            else
                load([resultPath,SegMethod,'MatrixZ/',fname],'Z');
            end
            lambdaLC = 1/lambda_lsr;
    end
    if plusLC
        fname = sprintf('%s_recovered_matrix_%d.mat',SegMethod,i);
        load([resultPath,SegMethod,'MatrixZ/',fname],'Z');
        Z = ssgd(X, K, lambdaLC, 4*K, Z, SegMethod, ProjMethod,0,0);
    end
    
    %% post processing
    if postProc
        disp('post processing ... ...');
        switch SegMethod
            case 'lrr'
                Z = rpcapsd(Z,0.8);

                % refining Z
                [U,S,V] = svd(Z);
                S = diag(S);
                r = min(4*K+1,sum(S>1e-3*S(1)));
                S = S(1:r);
                U = U(:,1:r)*diag(sqrt(S));
                U = normr(U);
                Z = U*U';Z=abs(Z);
                L = Z.^4.5;
                %         L = Z.^2;

                % spectral clustering
                L = (L + L')/2;
                D = diag(1./sqrt(sum(L,2)));
                L = D*L*D;
                CKSym = L;
                
            case 'ssc'
                rho = 0.7;
                C = 0.5*(abs(Z)+abs(Z'));
                CKSym = BuildAdjacency(thrC(C,rho));
        end
    else
        CKSym = 0.5*(abs(Z)+abs(Z'));
    end
    
    %% spectral clustering    
    grps = SpectralClustering(CKSym,K);
    err0 = Misclassification(grps,gnd);
    
    % show results
    disp(['seq ' num2str(i) ',err=' num2str(err0)]);
    errs(i) = err0;

end

%% show results
disp('results of all 156 sequences:');
disp(['max = ' num2str(max(errs)) ',min=' num2str(min(errs)) ...
    ',median=' num2str(median(errs)) ',mean=' num2str(mean(errs)) ',std=' num2str(std(errs))] );

errs = errs([1:id-1,id+1:end]);
disp('results of all 155 sequences:');
disp(['max = ' num2str(max(errs)) ',min=' num2str(min(errs)) ...
    ',median=' num2str(median(errs)) ',mean=' num2str(mean(errs)) ',std=' num2str(std(errs))] );

%% save results
fname = sprintf('errs_%s_%s_plus_%d_postproc_%d.mat',...
    SegMethod, ProjMethod,plusLC,postProc);
save([resultPath,fname],'errs');