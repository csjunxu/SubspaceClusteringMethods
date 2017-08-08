clc, clear all;

addpath(genpath(fullfile(pwd)));
cd '../../data/berkeley2M200';
% berkeley2M200

% trajectories associated with each motion live in
% an affine subspace of dimension at most three
% or a linear subspace of dimension at most four containing the affine
% subspace
p=1;
for trial = 1:1:100
    maxNumGroup = 5;
    for i = 1:maxNumGroup
        num(i) = 0;
    end
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
            cd ..
            
            if (foundValidData)
                n = max(s);
                N = size(x,2);
                F = size(x,3);
                D = 2*F; % 2F
                X = reshape(permute(x(1:2,:,:),[1 3 2]),D,N);
                
                if p == 0
                    X = X;
                else
                    p = 4*n;
                    X = DataProjection(X,p);
                end
                
                X = X';
                
                k = 3;    % scc 3;lssc 4, dimension of subspaces
                %[sampleLabels, averageL2Error] = lscc(X,k,n);
                [sampleLabels, averageL2Error] = scc(X,k,n);
                
                
                % compute misclassification rate
                SCCmissrate=missclass(sampleLabels,s,n)/length(s);
                disp(['Missclassification error: ' num2str(SCCmissrate*100) '%']);% averageL2Error
                
                num(n) = num(n) + 1;
                missrateTot{trial,n}(num(n)) = SCCmissrate;
                
                eval(['cd ' filepath]);
                save SCC_H.mat SCCmissrate
                cd ..
            end
        end
    end
    
    
    L = [2];
    for i = 1:length(L)
        j = L(i);
        avgmissrate(trial,j) = mean(missrateTot{trial,j});
        medmissrate(trial,j) = median(missrateTot{trial,j});
    end
end
for i = 1:length(L)
    j = L(i);
    allavgmissrate(j) = mean(avgmissrate(:,j));
    allmedmissrate(j) = median(medmissrate(:,j));
end
save SCC_B_3-4n.mat missrateTot avgmissrate medmissrate allavgmissrate allmedmissrate
