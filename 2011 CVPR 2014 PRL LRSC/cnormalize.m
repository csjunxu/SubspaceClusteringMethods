function [xn,normx] = cnormalize(x)
[K,N] = size(x);
normx = sqrt(sum(conj(x).*x,1));
xn =  x ./ (ones(K,1)*normx);
