clc;    clear;  close all;  warning off

addpath codes2
addpath codes

filename = 'YaleBCrop025_NVR3'; load(strcat(filename,'.mat'));
fid = fopen([filename '.txt'],'a+');

r = 6;
nSet = [2 3 5 8 10];
toler = 1e-3;   maxiter = 200;

lambda = 5; gamma1 = 1; eta1 = 0.0002;  gamma2 = 2; eta2 = 0.005;

for i = [ 1 : 5 ]
    
    n = nSet(i);    idx = Ind{n};   gnd = s{n};
    
    totacc = 0;     avgacc = 0;     Acc = zeros(size(idx,1),1);
    for j = 1 : size(idx,1)
        
        A = []; ph = 0;
        for p = 1 : n
            for h = 1:size(YY,3)
                ph = ph+1;  A(:,:,ph) = YY(:,:,h,idx(j,p));
            end
        end
        A = mat2gray(A);
        
        [Z,P,Q] = NVR3(A,r,lambda,gamma1,eta1,gamma2,eta2,maxiter,toler);
        
        acc = KSC_Acc(Z,2,n,gnd);
        
        Acc(j) = acc; totacc = totacc+acc;  avgacc = totacc/j;
        
        fprintf(1,'#Obj = %3d(%3d/%3d), acc: %6.4f, avgacc: %6.4f\n',n,j,size(idx,1),acc,avgacc);
        fprintf(fid,'#Obj = %3d(%3d/%3d), acc: %6.4f, avgacc: %6.4f\r\n',n,j,size(idx,1),acc,avgacc);
        
    end
    
end

fclose all;