clear;

% %% YaleB = 'YaleB.mat';
% load 'YaleBa.mat';

%% USPS
load 'C:\Users\csjunxu\Desktop\SC\Datasets\USPS.mat';
fea = fea' - min(min(fea));

% %% MNIST
% load 'MNIST.mat';
% fea = fea';



nSet = max(gnd);

Y = cell(nSet, 1);
S = cell(nSet, 1);
for i = 1:nSet
    Y{i, 1} = fea(:, gnd==i);
    S{i, 1} = gnd(find(gnd==i));
end
Ind = USPSInd;

save C:\Users\csjunxu\Desktop\SC\Datasets\USPS_Crop.mat Y S Ind;
    