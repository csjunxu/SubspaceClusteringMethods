clear;

Original_image_dir = 'C:\Users\csjunxu\Desktop\SC\Results\';
fpath = fullfile(Original_image_dir, 'YaleB_LSR*.mat');
mat_dir  = dir(fpath);
mat_num = length(mat_dir);

meanavgmissrate = [];
for i = 1 : mat_num
    fprintf([mat_dir(i).name '\n']);
    eval(['load C:\Users\csjunxu\Desktop\SC\Results\' mat_dir(i).name]);
    meanavgmissrate = [meanavgmissrate mean(avgmissrate)];
end

[minmeanavgmissrate, index] = sort(meanavgmissrate, 'ascend');
fprintf([num2str(minmeanavgmissrate(1)) '\n']);
fprintf('The minimum one is %s', mat_dir(index(1)).name);
