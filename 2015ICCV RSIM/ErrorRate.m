function Missrate = ErrorRate(Seg,RefSeg)

if(max(max(RefSeg))~=1 && max(max(Seg))~=1) %Calculate the error rate using the results of standard spectral clustering method
	%Missrate = Misclassification(Seg,RefSeg);
	Missrate = missclassGroups( Seg,RefSeg,max(RefSeg) ) ./ length(RefSeg);
elseif(max(max(RefSeg))==1 && max(max(Seg))==1)
	[N,n] = size(RefSeg); %number of clusters
	SegVec = zeros(N,1);
	RefSegVec = zeros(N,1);
	for i=1:n
		SegVec = SegVec+i*Seg(:,i);
		RefSegVec = RefSegVec+i*RefSeg(:,i);
	end
	%Missrate = Misclassification(SegVec,RefSegVec);
	Missrate = missclassGroups( SegVec,RefSegVec,max(RefSegVec) ) ./ length(RefSegVec);
elseif(max(max(RefSeg))==1 && max(max(Seg))~=1)
	[N,n] = size(RefSeg); %number of clusters
	RefSegVec = zeros(N,1);
	for i=1:n
		RefSegVec = RefSegVec+i*RefSeg(:,i);		
	end
	%Missrate = Misclassification(Seg,RefSegVec);
	Missrate = missclassGroups( Seg,RefSegVec,max(RefSegVec) ) ./ length(RefSegVec);
else
	[N,n] = size(Seg); %number of clusters
	SegVec = zeros(N,1);
	for i=1:n
		SegVec = SegVec+i*Seg(:,i);		
	end
	%Missrate = Misclassification(SegVec,RefSeg);	
	Missrate = missclassGroups( SegVec,RefSeg,max(RefSeg) ) ./ length(RefSeg);
end
