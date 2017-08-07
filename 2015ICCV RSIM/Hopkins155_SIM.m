% Motion segmentation with Robust Shape Interaction Matrix Method
% Pan Ji, pan.ji@anu.edu.au
% Nov 2014, @ANU
addpath(genpath(pwd))
cd('C:/Users/csjunxu/Desktop/SC/Datasets/Hopkins155/') % cd to your Hopkins155 file path
clear; close all

addpath('C:\Users\csjunxu\Desktop\SC\Ncut_9');

warning off

file = dir;
ii = 0;
ii2 = 0;
ii3 = 0;
tic
for i = 1:length(file)
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
				if(max(s)==5)
					foundValidData = false;
				end
                break
            end
		end        
        cd ..
		
		if(foundValidData)
			N = size(x,2);
			F = size(x,3);
			D = 3*F;						
									
			X = reshape(permute(x(1:3,:,:),[1 3 2]),D,N);	% note here the all-one rows are also included
						
			[missrate, grp, bestRank, minNcutValue,W] = RSIM(X, s);			

			ii = ii+1;	
			Missrate(ii) = missrate;					
			disp([filepath ': ' num2str(100*Missrate(ii)) '%, dim:' num2str(bestRank) ', nMotions: ' num2str(max(s)) ', seq: ' num2str(ii)]);
			if(max(s)==2)
				ii2 = ii2+1;
				Missrate2(ii2) = Missrate(ii);				
			else
				ii3 = ii3+1;
				Missrate3(ii3) = Missrate(ii);				
			end
		end
	end
end
time = toc;
avgtime = time/ii

avgtol = mean(Missrate);
medtol = median(Missrate);
avgtwo = mean(Missrate2);
medtwo = median(Missrate2);
avgthree = mean(Missrate3);
medthree = median(Missrate3);

disp('Results on Hopkins155')
disp(['Mean of all: ' num2str(100*avgtol) '%' ', median of all: ' num2str(100*medtol) '%;']);
disp(['Mean of two: ' num2str(100*avgtwo) '%' ', median of two: ' num2str(100*medtwo) '%;']);
disp(['Mean of three: ' num2str(100*avgthree) '%' ', median of three: ' num2str(100*medthree) '%.']);