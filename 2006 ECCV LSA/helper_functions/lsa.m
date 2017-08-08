function [group,ranks]=lsa(x,n,k,d,constmodelselection,spectralalg)
%function [group,ranks]=lsa(x,n,k,d,constmodelselection,spectralalg)
%
%   Local Subspace Affinity
%       x                   data points
%       n                   number of planes
%       d                   dimension of the planes (default, dimension of x minus one)
%       k                   neighbours to use (default or if k==0, k=d)
%       constmodelselection if set to >0, use model selection techniques to
%                           find the dimension of the subspaces (default 0)
%       spectralalg         select the spectral clustering algorithm to use:
%           'kmeans'     Ng's algorithm (default)
%           'normalcuts' Normalized cuts   

%   Ver. 4/5/2012

DEBUG=0;

if(exist('constmodelselection','var')==0)
    constmodelselection=0;
end

if(exist('spectralalg','var')==0)
    spectralalg='kmeans';
end

[K,N]=size(x);

if(exist('d','var')==0)
    d=K-1;
end

if(exist('k','var')==0 || k==0)
    k=d;
end


x = cnormalize(x);

%find neighbours
angles=acos(min(max(x'*x,-1),1));
[sorted,neighbours]=sort(angles);

%estimate local normal
bases=zeros(size(x,1),d,N);
ranks=zeros(N,1);
for(i=1:N)
    [U,S,V]=svd(x(:,neighbours(1:k+1,i)));
    if(constmodelselection>0)
        ranks(i)=max(modelselection(diag(S),constmodelselection),2);
        bases(:,1:ranks(i),i)=U(:,1:ranks(i));
        if(DEBUG==1)
            disp(['Rank point ' num2str(i) ':' num2str(ranks(i))])
        end
    else
        bases(:,:,i)=U(:,1:d);
        ranks(i)=d;
    end
end

%compute angle between subspaces => similarity matrix
for(j=1:N)
     for(i=j:N)
         if(constmodelselection>0)
             distance(j,i) = subspace(bases(:,1:ranks(j),j),bases(:,1:ranks(i),i));
             %distance(j,i) = subspaceaffinity(bases(:,1:ranks(j),j),bases(:,1:ranks(i),i));
         else
             distance(j,i) = subspace(bases(:,:,j),bases(:,:,i));
             %distance(j,i) = subspaceaffinity(bases(:,:,j),bases(:,:,i));
         end
         distance(i,j) = distance(j,i);
     end
end
if(DEBUG==1)
    imshow(abs(distance),[]);
end

%spectral clustering
switch(spectralalg)
    case 'kmeans'
        [diagMat,LMat,X,Y,group,errorsum]=spectralcluster(distance,n,n);
    case 'normalcuts'
        group=spectralclusternormalcut_recursive(n,distance);
end
