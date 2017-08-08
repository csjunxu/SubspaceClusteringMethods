% Making the diagonal line of coefficient matrix equals to zeros.
function [Ar] = avoidze(datam)
[nda, mda] = size(datam);
for j = 1 : mda
    Da = datam;
    datam(:,j) = 0;
    Ar(:,:,j) = datam;
    datam = Da;
end