%% THE ALC ALGORITHM FOR MOTION SEGMENTATION
%
%
%
% The website is 'http://perception.csl.illinois.edu/coding/motion/'
% Modified by Xujun

%%
clc, clear all

addpath 'C:/Users/csl/Desktop/motionseg/CodeALC/helpers';
addpath 'C:/Users/csl/Desktop/motionseg/CodeALC';
sequence_Dir = 'C:/Users/csl/Desktop/motionseg/Hopkins155';
cd 'C:/Users/csl/Desktop/motionseg/Hopkins155';

maxNumGroup = 5;
for i = 1:maxNumGroup
    num(i) = 0;
end

epsilon = logspace(-5,3,101);

d = dir;
for i = 1:length(d)
    if ( (d(i).isdir == 1) && ~strcmp(d(i).name,'.') && ~strcmp(d(i).name,'..') )
        filepath = d(i).name;
        eval(['cd ' filepath]);
        
        f = dir;
        foundValidData = false;
        for j = 1:length(f)
            if ( ~isempty(strfind(f(j).name,'_truth.mat')) )
                ind = j;
                foundValidData = true;
                break
            end
        end
        eval(['load ' f(ind).name]);
        n = max(s);
        cd ..
        if (foundValidData)            
            sequence_Name = [filepath '/' f(ind).name];
            [rawData, trueLabels] = load_sequence(sequence_Dir,sequence_Name);
            processedData = process_sequence(rawData, true);
            result = try_sequence(sequence_Name, processedData, epsilon);
            computedLabels = find_best_segmentation(result, processedData, n, epsilon);
            err = compare_labels(trueLabels, computedLabels);
            
            num(n) = num(n) + 1;
            errTot{n}(num(n)) = err;
            
            eval(['cd ' filepath]);
            save ALC_MS.mat err result
            cd ..
        end
    end
end

L = [2 3];
for i = 1:length(L)
    j = L(i);
    avgmissrate(j) = mean(errTot{j});
    medmissrate(j) = median(errTot{j});
end

save ALC_MS.mat errTot avgmissrate medmissrate