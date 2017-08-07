1. Download Hopkins motion segmentation dataset from http://vision.jhu.edu/data/.

2. Install Ncuts function from https://www.cis.upenn.edu/~jshi/software/ .

% Obviously the new_eigs() function in ncuts is incompatible with the ARPACK 
% version in the latest Matlab
% Solution£º
% Use Matlab's eigs() function instead of the eigs_new() provided in the normalized 
% cuts package. I guess eigs_new() was designed to solve some compatibility issue with 
% a previous version of Matlab, and is now itself causing an issue.

3. Change filepath where applicable.