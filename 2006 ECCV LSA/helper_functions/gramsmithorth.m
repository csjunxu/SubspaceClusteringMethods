%function y=gramsmithorth(x)
%   Returns Y the Gram-Smith orthogonalization of the colums of X
%   Y and X have the same dimensions
function y=gramsmithorth(x)
[K,D]=size(x);

I=eye(K);

y=x(:,1)/norm(x(:,1));
for(i=2:D)
    newcol=(I-y*y')*x(:,i);
    newcol=newcol/norm(newcol);
    y=[y newcol];
end
