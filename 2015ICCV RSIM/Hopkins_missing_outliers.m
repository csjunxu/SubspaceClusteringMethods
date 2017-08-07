% Pan Ji, pan.ji@anu.edu.au
% Nov 2014, @ANU
addpath(genpath(pwd))
cd('D:\Reference_Code\Hopkins155_AdditionalSequences_OutliersAndMissingData') 
% Download the dataset from http://vision.jhu.edu/data/
clear; close all

file = dir;
ii = 0;
for i = 9:length(file)
	if( (file(i).isdir == 1) && ~strcmp(file(i).name,'.') && ~strcmp(file(i).name,'..') )
		filepath = file(i).name;
		eval(['cd ' filepath]);
		
		f = dir;
        foundValidData = false;
        for j = 1:length(f)
            if( ~isempty(strfind(f(j).name,'_truth.mat')) )
                ind = j;
                foundValidData = true;
				eval(['load ' f(ind).name]);
                break
            end
		end        
        cd ..
		
		if(foundValidData)
			N = size(x,2);
			F = size(x,3);
			D = 2*F;	
						
			X = reshape(permute(x(1:2,:,:),[1 3 2]),D,N);
			
			indX = zeros(N,D);
			for jj = 1:F				
				indX(:,2*jj-1) = m(:,jj);
				indX(:,2*jj) = m(:,jj);
			end
			idx_all0 = [];
			for kk=1:2*F
				if(isequal(indX(kk,:),zeros(1,2*F)))
					idx_all0 = [idx_all0,kk];
				end
			end
			indX = indX';
			X(:,idx_all0) =[];
			indX(:,idx_all0) =[];
			s(idx_all0) = [];			
			
			X = X.*indX;			
			
			[missrate, grp, bestRank, minNcutValue] = RSIM_Incomplete(X',indX',s,4,1);

			ii = ii+1;
			Missrate(ii) = missrate;
			disp([filepath ': ' num2str(100*Missrate(ii)) '%' ', nMotions: ' num2str(max(s)) ', seq: ' num2str(ii)]);
		end
	end
end

avgtol = mean(Missrate);
medtol = median(Missrate);
maxtol = max(Missrate);
stdtol = std(Missrate);

disp('Results on Hopkins Missing data')
disp(['Mean: ' num2str(100*avgtol) '%' ', median: ' num2str(100*medtol) '%;'...
      ', max: ' num2str(100*maxtol) '%;' ', std: ' num2str(100*stdtol) '%;']);
