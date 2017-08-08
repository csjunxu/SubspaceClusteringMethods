% find the nth column coefficient of matrix
function[finalXn] = finalresulXn(frxdata,K,firwei,Sig,lam,n,X)
       xnl = size(frxdata,2);
      idlam = 2*lam*eye(xnl);
      [fxArgdat] = avoidze(frxdata);
      [fpartXn] = solvematrix(frxdata,K,firwei,Sig,n,X);
      finalXn = inv(fpartXn*fxArgdat(:,:,n) + idlam)*(fpartXn*frxdata(:,n));
     end