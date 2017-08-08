function acc = KSC_Acc(X,aa,K,gnd)

[U, S, V]=svd(X);
S = diag(S);
r = sum(S>1e-6);
S = diag(S(1 : r));
U = U(:, 1 : r);
V = V(:, 1 : r);

M = U * S.^(1/2);
% M = U * diag(diag(S).^(1/2));
mm = normr(M);
rs = mm * mm';
for t=1:length(aa)
    L = rs.^(2 * aa(t));
    Idx = spectral_clustering(L, K);
    %     acc = 1 - calculate_err(gnd, actual_ids);
%     acc = 1 - ClusterErr(K,gnd,Idx);
    acc = accuracy(gnd,Idx)/100;
end

end