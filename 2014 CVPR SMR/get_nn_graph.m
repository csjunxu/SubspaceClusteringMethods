function [pairs,wcost,numpairs]=get_nn_graph(Xp,knn)

[~,P] = size(Xp);
neighborPoints =knn;

distMin = squareform(pdist(Xp'));


[minT,midx] = sort(distMin,'ascend');
neighMat = zeros(P,P);
for pp = 1:P
    neighMat(midx(1:neighborPoints,pp),pp) = 1;
    neighMat(pp,midx(1:neighborPoints,pp)) = 1;
end

neighMat = tril(neighMat,-1);

[nzidx] = find(neighMat>0);
[nzr,nzc] = ind2sub(size(neighMat),nzidx);
pairs = [nzr nzc]'-1;

numpairs = length(nzidx);

wcost = ones(length(nzidx),1);
