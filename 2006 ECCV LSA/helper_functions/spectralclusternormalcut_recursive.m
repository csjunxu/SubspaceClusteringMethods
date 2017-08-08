%function group=spectralclusternormalcut_recursive(n,simMat)
%n      final number of cluster
%simMat similarity matrix
function group=spectralclusternormalcut_recursive(n,simMat)
%trivial case with 1 cluster
if(n==1)
    group=ones(size(simMat,1),1);
end
%initial bipartition with spectral clustering
group=spectralclusternormalcut(simMat);

%call spectral clustering recursively
for(i=2:n-1)
    mincheeger=Inf;     %init
    %try to split each cluster
    for(j=1:i)
        clusterindeces=find(group==j);
        [groupcluster,clustercheeger]=spectralclusternormalcut(simMat(clusterindeces,clusterindeces));
        %if the cheeger constant is less than in the previous clusters
        if(clustercheeger<mincheeger)
            %save the new partition
            mincheeger=clustercheeger;
            minclusterindeces=clusterindeces;
            mingroupcluster=groupcluster;
        end
    end
    %update with the new partition
    group(minclusterindeces(find(mingroupcluster==2)))=i+1;
end
