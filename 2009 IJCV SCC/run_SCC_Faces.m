clear, close,clc;

load 'YaleBCrop025.mat';

nSet = [2 3 4 5 6 7 8 9 10];
for i = 1:length(nSet)
    n = nSet(i);
    idx = Ind{n};
    disp(['Clustering ' num2str(n) ' faces!']);
    for trial = 1:1:100
        disp([' The ' num2str(trial) 'th trial!']);
        for j = 1: size(idx,1)
            X = [];
            for p = 1:n
                X = [X Y(:,:,idx(j,p))];
            end
            [D,N] = size(X);
            
            X = X';
            
            k = 3;    % scc 3;lssc 4, dimension of subspaces
            [sampleLabels, averageL2Error] = scc(X,k,n);
            %[sampleLabels, averageL2Error] = lscc(X,k,n);
            
            
            % compute misclassification rate
            missrate=missclass(sampleLabels,s{n},n)/length(s{n});
            disp(['Missclassification error: ' num2str(missrate*100) '%']);% averageL2Error
            
            missrateTot{trial,n}(j) = missrate;
        end
        avgmissrate(trial,n) = mean(missrateTot{trial,n});
        medmissrate(trial,n) = median(missrateTot{trial,n});
    end
    allavgmissrate(n) = mean(avgmissrate(:,n));
    allmedmissrate(n) = median(medmissrate(:,n));
    save SCC_Faces.mat missrateTot avgmissrate medmissrate
end
save SCC_Faces.mat missrateTot avgmissrate medmissrate allavgmissrate allmedmissrate