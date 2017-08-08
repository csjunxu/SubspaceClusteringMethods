% file svd_thin_US_rot.m
% (c) 2008 Stefan Atev
% A faster, non-truncating implementation of M. Brand's incremental SVD
% The subspace update is novel, and note that Ur is not
% orthogonal.
function [U, S]= svd_thin_US_rot(X, r, tol)
    [p, q]= size(X); % p rows, q columns
    er= 0; % effective rank is 0 before any columns are appended
    U= [];
    Ur= [];
    invUr= [];
    S= [];
    if nargin< 3
        tol= 1e-10;
    end
    for i= 1:q
        a= X(:, i);
        if er== 0
            s= norm(a);
            if s> tol
                U= a/ s;
                Ur= 1;
                invUr= 1;
                S= s;
                er= 1;
            end
        else
            m= Ur'* (U'* a);
            p= a- U* (Ur* m);
            np= norm(p);
            if (np> tol)
                % non-negligible residual
                if (er< r)
                    % add a column to U
                    M= [S, m; zeros(1, er), np];
		      if (sum(isnan(M(:))) > 0)
                       M(isnan(M)) = 0;
                    end
                    [mU, S]= svd(M);
                    AB= mU(1:er,1:er+1);
                    CD= mU(er+1,1:er+1);
                    Ur= [Ur* AB; CD];
                    invUr= [AB'* invUr, CD'];
                    U= [U, p/ np];
                    er= er+ 1;
                else
                    % cannot add column, attempt to update U
                    % in O(er^2) time
                    M= [S, m; zeros(1, er), np];
		      if (sum(isnan(M(:))) > 0)
                        M(isnan(M)) = 0;
                    end
                    [mU, S]= svd(M);
                    A= mU(1:er,1:er);
                    B= mU(1:er,er+1);
                    C= mU(er+1,1:er);
                    D= mU(er+1,er+1);
                    if abs(D)< eps
                        % fast update will fail, do O(er^3) update
                        AB= mU(1:er,1:er+1);
                        CD= mU(er+1,1:er+1);
                        U= U* (Ur* AB)+  (p* CD/ np);
                        U= U(:,1:er);
                        Ur= eye(er);
                        invUr= Ur;
                    else
                        % fast update can work
                        invA= (A- B* C/ D)';
                        Ur= Ur* A;
                        invUr= invA* invUr;
                        mm= C* invUr;
                        U= U+ p* (mm/ np);
                    end
                    S= S(1:er,1:er);
                end
            else
                % truncating update, residual is ignored
                M= [S, m];
		  if (sum(isnan(M(:))) > 0)
                   M(isnan(M)) = 0;
                end
                [mU, mS]= svd(M);
                Ur= Ur* mU(1:er,1:er);
                invUr= mU(1:er,1:er)'* invUr;
                S= mS(1:er,1:er);
            end
        end
    end
    U= U* Ur;
end
