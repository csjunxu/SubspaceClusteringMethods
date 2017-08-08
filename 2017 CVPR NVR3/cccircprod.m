function Y = cccircprod(X,p)
[a,b,n] = size(X);
[~,r] = size(p);
Y = zeros(a*r,n);
for t = 1:n
    y = X(:,:,t)*p;
    Y(:,t) = reshape(y,a*r,1);
    %     Y(:,t) = X(:,:,t)*p;
end

end