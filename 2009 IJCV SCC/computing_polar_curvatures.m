function curvatures = computing_polar_curvatures(S,sampledColumns,d)
% The strategy used in this code for fast computation of polar curvatures 
% is based on Stefan Atev's implementation of a kernelized version.
% However, we have improved the special case (kernel = regular dot product),
% which is what we need here.

m = size(S,1);
c = size(sampledColumns,2);

% build an index of which Gramian columns will be needed
%sampledPoints= union(reshape(sampledColumns,1,[]),[]);
sampledPoints= unique(sampledColumns).';
numberSampledPoints= numel(sampledPoints);

colMap= zeros(1, m);
colMap(sampledPoints)= 1:numberSampledPoints;

K = S*S(sampledPoints,:)'; % m*c dot product matrix

norms = zeros(m,1);
norms(sampledPoints) = diag(K(sampledPoints,:));

others = setdiff(1:m,sampledPoints);
norms(others) = sum(S(others,:).^2,2);

curvatures = Inf(m,c);
%curvatures = zeros(m,c);

if d > 0
    
    for col= 1:c
        % computes one full column of polar curvatures directly 

        Fidx= sampledColumns(:,col); % Fidx = (j_1,...,j_d+1)
        Ridx= setdiff(1:m, Fidx); % R_idx = 1:m - (j_1,...,j_d+1)

        K_FF= K(Fidx, colMap(Fidx));
        K_RF= K(Ridx, colMap(Fidx));
        k_RR= norms(Ridx);

        % compute volume for each Ridx
        [V, L]= eig(K_FF+ ones(d+1));
        L= diag(L);
        detK_FF= prod(L);
        adjK_FF= V* diag(detK_FF./ L)* V';
        K_RF1= K_RF+ 1;
        volumes= (k_RR+ 1)* detK_FF- sum((K_RF1* adjK_FF).* K_RF1, 2);

        % compute edge lengths for each Ridx
        B_FF= repmat(norms(Fidx), [1 (d+ 1)])+ repmat(norms(Fidx)', [(d+ 1) 1])- 2*K_FF;
        B_RF= repmat(k_RR, [1 (d+ 1)])+ repmat(norms(Fidx)', [(m- d- 1) 1])- 2* K_RF;
        diams = max([B_RF repmat(max(max(B_FF)),m-d-1,1)],[],2);        
        
        plen= [repmat(prod(eye(d+ 1)+ B_FF), [(m- d- 1) 1]).* B_RF, prod(B_RF, 2)];
        curvatures(Ridx,col)= diams.*sum(repmat(volumes,1,d+2)./plen,2);
        
    end

else
    % special case d= 0, just needs edge norms
    for col= 1:c
        
        Fidx= sampledColumns(1,col);
        Ridx= setdiff(1:m, Fidx);
        K_RF= K(Ridx, colMap(Fidx));
        k_RR= norms(Ridx);
        curvatures(Ridx,col)= (k_RR+ norms(Fidx))- 2*K_RF;
        
    end
    
end

%curvatures = sqrt(abs(curvatures));
curvatures = abs(curvatures);
curvatures(isnan(curvatures))=0;
