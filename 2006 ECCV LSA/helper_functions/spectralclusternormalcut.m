function [group,cheeger]=spectralclusternormalcut(simMat)
%number of thresholds to try
l=1500;

d=sum(simMat,2);
D=diag(d);
D2=diag(sqrt(d));
[U,S,V]=svd(D2*(D-simMat)*D2);
x=U(:,end-1);
th=linspace(min(x),max(x),l+2);
th=th(2:end-1);
cutcosts=zeros(l,1);
cheegers=zeros(l,1);
for(i=1:l)
    cutcosts(i)=evaluatenormalcut(x>th(i),simMat);
end
[mincost,minindex]=min(cutcosts);
group=(x>th(minindex))+1;
cheeger=cheegerpartition(group,simMat);
