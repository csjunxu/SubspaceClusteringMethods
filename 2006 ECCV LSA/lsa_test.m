clc, clear all

addpath(genpath(fullfile(pwd,'/helper_functions')));
cd '../../data/berkeley2M200';

ds = 'B';

if ds =='H'
        L = [2 3];
    elseif ds=='B'
        L = [2];
    end
maxNumGroup = 5;
for i = 1:maxNumGroup
    num(i) = 0;
end

for k = 3  %dimension of subspaces: 4 for rigid while 7 for non-rigid
    for kNeigh= [30];  %number of neighbors
        for  p= [0]
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
                        D = 2*F;
                        X = reshape(permute(x(1:2,:,:),[1 3 2]),D,N);
                        
                        % projection
                        if p == 0
                            cas = '2F';
                        else
                            p = 4*n;
                            cas = '4n';
                        end
                        X = DataProjection(X,p);
                        
                        %call LSA
                        group=lsa(X,n,kNeigh,k);
                        
                        %compute misclassification rate
                        LSAmissrate=missclassGroups(group,s,n)/length(s);
                        fprintf('seq %s\t %f\n', d(i).name , LSAmissrate);
                        
                        num(n) = num(n) + 1;
                        missrateTot{n}(num(n)) = LSAmissrate;
                        
                        eval(['cd ' filepath]);
                        save LSA_MS.mat LSAmissrate
                        cd ..
                    end
                end
            end
            for i = 1:length(L)
                j = L(i);
                avgmissrate(j) = mean(missrateTot{j});
                medmissrate(j) = median(missrateTot{j});
            end
            name = sprintf('LSA_%s_%s_k%dkNeigh%d.mat',cas,ds,k,kNeigh);
            save(name, 'missrateTot','avgmissrate','medmissrate','k','kNeigh');
        end
    end
end
