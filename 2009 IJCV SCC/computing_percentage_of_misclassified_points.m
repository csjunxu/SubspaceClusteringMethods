function p = computing_percentage_of_misclassified_points(indices,trueLabels)

[sortedLabels, inds_sort] = sort(trueLabels, 'ascend');
indices = indices(inds_sort);

N = numel(trueLabels);
K = sortedLabels(end);

%if K == 1; p = 0; return;

planeSizes = zeros(K,1);
i = 1; ini = 1;
for k = 1:K-1,
    while sortedLabels(i) == k; i = i+1; end;
    planeSizes(k) = i - ini;
    ini = i;
end
planeSizes(K) = N+1 - ini;

num = zeros(K,K);
for k = 1:K
    for j = 1:K
        num(k,j) = sum(indices(sum(planeSizes(1:k-1))+1:sum(planeSizes(1:k)))==j);
    end
end

p = 1-maximum_number_of_correctly_classified_points(num)/sum(planeSizes);

%%

function n = maximum_number_of_correctly_classified_points(num)

K = size(num,1);

if K > 2 
    n = zeros(K,1);
    for j = 1:K
        n(j) = num(1,j)+maximum_number_of_correctly_classified_points(num(2:end,[1:j-1 j+1:K]));
    end
    n = max(n);
else
    n = max(num(1,1)+num(2,2), num(1,2)+num(2,1));
end