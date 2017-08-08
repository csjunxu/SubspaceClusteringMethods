warning off all
data_dir = './data/usps';
addpath(genpath(data_dir));
load usps_part.mat
exp_name = 'HandwrittenDigit';
data = fea';

para.knn = 4;
para.gamma = 5;
para.elpson = 0.001;
para.alpha = 2.^-16;

para.aff_type = 'J1';

W = smr(data,para);
W2 = W;
for ic = 1 : size(W,2)
   W2(:,ic) = W(:,ic) / (max(abs(W(:,ic)))+eps) ;    
end
[groups] = clu_ncut(W2,max(gnd));
ce = compacc_ce(groups,gnd)

