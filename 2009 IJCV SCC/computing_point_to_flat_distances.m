function distancePoints2Flats = computing_point_to_flat_distances(X,centers,bases)

[N,D] = size(X);
K = length(centers);

distancePoints2Flats = Inf(N, K);
for k = 1:K

    if ~isempty(centers{k}) && ~isempty(bases{k})
        
        X_centered = X - repmat(centers{k},N,1);
        distancePoints2Flats(:,k) = ...
            sum(X_centered.^2,2) - sum((X_centered*bases{k}').^2,2);
            %sum(((X - repmat(centers{k},N,1))*(eye(D) - bases{k}'*bases{k})).^2,2);
    end
    
end