%evaluates the normal cut function
%   group      is a vector of zeros and ones that indicates the two partitions
%   simMat     is the similarity matrix
function cost=evaluatenormalcut(group,simMat);

d=sum(simMat,2);                        %grade of each node (sum of distances on the row)
assocA=sum(d(find(group==0)));          %association of cluster A
assocB=sum(d(find(group==1)));          %association of cluster B
[IcutA,IcutB]=meshgrid(group,1-group);  %bool that indicates if a group is connected to A and/or B
IcutAB=and(IcutA,IcutB);                %indeces of nodes of A connected to B
IcutBA=not(or(IcutA,IcutB));            %indeces of nodes of B connected to A
cutAB=sum(simMat(find(IcutAB)));        %sum of arcs from A to B
cutBA=sum(simMat(find(IcutBA)));        %sum of arcs from B to A
cost=cutAB/assocA+cutBA/assocB;         %normal cut cost function
