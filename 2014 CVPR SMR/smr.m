function [L] = smr(X,para)
%This routine solves the smooth representation problem with fro-norm data term
%
alpha = para.alpha;
gamma = para.gamma;

[d,n] =size(X);
[pairs,wcost,numpairs]=get_nn_graph(X,para.knn);

nX = sqrt(sum(X.^2));

R = zeros(n,numpairs);
for i=1:numpairs
    R(pairs(1,i)+1,i) = wcost(i);
    R(pairs(2,i)+1,i) = -wcost(i);
end

R = R/(para.knn-1);

xtx = X'*X;
rtr = 0.5*R*R';
% elpson = 0.01;
A = alpha*xtx;
B = rtr+para.elpson*eye(size(xtx));
C = -alpha*xtx;
J = lyap(A,B,C); 
if strcmp(para.aff_type,'J1')
    L =(abs(J)+abs(J'))/2;
elseif strcmp(para.aff_type,'J2')
    L=abs(J'*J./(nX'*nX)).^gamma;
elseif strcmp(para.aff_type,'J2_nonorm')
    L=abs(J'*J).^gamma;
end


