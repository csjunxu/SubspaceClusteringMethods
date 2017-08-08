%evaluates the cheeger constant for a given partition
function h=cheegerpartition(group,simMat);
d=sum(simMat,2);                        %grade of each node (sum of distances on the row)
[IcutA,IcutB]=meshgrid(group-1,2-group);  %bool that indicates if a group is connected to A and/or B
IcutAB=and(IcutA,IcutB);                %indeces of nodes of A connected to B
cutAB=sum(simMat(find(IcutAB)));        %sum of arcs from A to B

%compute the cheeger constant for this partition
h=cutAB/min(sum(d(find(group==1))),sum(d(find(group==2))));