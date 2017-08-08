close all
clear all
clc
tic
disp('Runing Do not close')
load('AR_Gray.mat');
load('groudtruth.mat');
contnu = 260;
new = AR_Gray(:,1:contnu);
s = groudtruth(1: contnu,:);
%--------------------------------------------
K = 5;
dedim= 31;
lamdmog = .8;% the shrinkage parameter
[eigvector, eigvalue, elapse] = PCA(new);
Dat = eigvector'*new;
minc = min(s);
maxc = max(s);
classs = maxc - minc +1;
%-------------------------------------------
Dat = Dat/norm(Dat);%normalize
Dat = Dat(1 :dedim,:);
W0 =1/K*ones(1,K);% initialize the weight
rorSi = size(Dat,1);
for j = 1 : K
    Sig0(:,:,j) =eye(rorSi);% initialize the covariance
end
flagmog = 1;
lop = size(Dat,2);
%---------------------------------
X0 = zeros(lop);%initialize the solution
%----------------------------------------------------------
jishuloop = 0;
err = 1e-4;
iterm = 1e+4;
while flagmog ==1
    jishuloop = 1 +jishuloop
    X1 = X0;
    XDat = Dat;
    for i = 1 : lop
        [finalXi] = finalresulXn(XDat,K,W0,Sig0,lamdmog,i,X0);
        %[xPro] = Mltg(xArgdat(:,:,n)*X(:,n)-xdata(:,n),Sig(:,:,j));
        X0(:,i) = finalXi;
    end
    Sig1 = Sig0;
    fcodat = Dat;
    for j = 1 : K
        [Sigmj] = findcovj(fcodat,Sig0,K,X0,W0,j);
        Sig0(:,:,j) = Sigmj ;
    end
    W1 = W0;
    wdatam = Dat;
    for j = 1 : K
        [newweij] = updattweig(wdatam,W0,j,K,X0,Sig0);
        W0(j) = newweij;
    end
    compnor = norm(X0 - X1)
    %-----------------------------------
    if compnor < err
        flagmog = 0
    else
        if jishuloop > iterm
            disp('Not converge')
            break
        end
    end
    %--------------------------------------
end
[idx] = clu_ncut(abs(X0) + abs(X0'),classs);
[ODX0,sc,OutlierIndx,Fail] = OutlierDetection(X0,s');
Ngu = max(size(Dat));
Adx0 = BuildAdjacency(ODX0,Ngu);%build the affine matrix
imshow(Adx0)
classs = double(classs);
[GAdx0, SingVals, LapKernel] = SpectralClustering(Adx0,classs);% classs:the number of subspaces
%------------------------------------------------
Missrate = Misclassification(GAdx0,sc);
accmog = max(1 - Missrate)
toc
disp('over')
