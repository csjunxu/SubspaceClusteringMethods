function error = evalSSR_error(C, s)
% C \in R^N-by-N: representation matrix by self-expressiveness based method
% s \in {1, 2, ... n}^N: group labels
%     
% error: average SSR error.

N = length(s);

error = 0;
e_vec = zeros(1, N);
for iN = 1:N
    if norm(C(:, iN), 1) < eps
        error = error + 1;
    else
        error = error + norm( C(s ~= s(iN), iN), 1 ) / norm(C(:, iN), 1);
    end
    e_vec(iN) = norm( C(s ~= s(iN), iN), 1 ) / norm(C(:, iN), 1);
end
error = error / N;

end
