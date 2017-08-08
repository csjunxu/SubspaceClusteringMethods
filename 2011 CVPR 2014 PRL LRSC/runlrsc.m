
clc, clear all

addpath(genpath(fullfile(pwd)));
cd '../../data/berkeley2M200';
% berkeley2M200
% Hopkins155
ds = 'B';

for tau  =232;
    
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
                D = 2*F;
                X = reshape(permute(x(1:2,:,:),[1 3 2]),D,N);
                
                r = 0; outlier = false; rho = 0.7;
                [missrate1,C1] = lrsc(X,tau,r,outlier,rho,s);
                disp(['Missclassification error :' num2str(missrate1*100) '%']);
                
                r = 4*n; affine = true; outlier = false; rho = 0.7;
                [missrate2,C2] = lrsc(X,tau,r,outlier,rho,s);
                disp(['Missclassification error :' num2str(missrate2*100) '%']);
                
                num(n) = num(n) + 1;
                missrateTot1{n}(num(n)) = missrate1;
                missrateTot2{n}(num(n)) = missrate2;
                
                eval(['cd ' filepath]);
                cd ..
            end
        end
    end
    
    if ds =='H'
        L = [2 3];
    elseif ds=='B'
        L = [2];
    end
    allmissrate1 = [];
    allmissrate2 = [];
    for i = 1:length(L)
        j = L(i);
        avgmissrate1(j) = mean(missrateTot1{j});
        medmissrate1(j) = median(missrateTot1{j});
        avgmissrate2(j) = mean(missrateTot2{j});
        medmissrate2(j) = median(missrateTot2{j});
        allmissrate1 = [allmissrate1 missrateTot1{j}];
        allmissrate2 = [allmissrate2 missrateTot2{j}];
    end
    avgallmissrate1 = sum(allmissrate1)/length(allmissrate1);
    medallmissrate1 = median(allmissrate1);
    avgallmissrate2 = sum(allmissrate2)/length(allmissrate2);
    medallmissrate2 = median(allmissrate2);
    
    name = sprintf('LRSC_%s_%d.mat',ds,tau);
    save(name,'missrateTot1','missrateTot2','avgmissrate1','medmissrate1','avgmissrate2','medmissrate2','avgallmissrate1','medallmissrate1','avgallmissrate2','medallmissrate2','tau');
end