% Calculating the Gaussian density function 
function [Pro] = Mltg(x,si)
    nsi = size(si,1);
    con = 1/sqrt(abs(det(si) + realmin)*(2*pi)^nsi);
    Pro = con*exp(-.5*x'*inv(si)*x);
end