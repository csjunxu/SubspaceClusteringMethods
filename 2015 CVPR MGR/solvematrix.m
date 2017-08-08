% Finding the nth column of coefficient matrix
function[partXn] = solvematrix(xdata,K,xwei,Sig,n,X)
       partXn = 0;
       [xArgdat] = avoidze(xdata);
        [xweisum] = weightsumgaupdf(Sig,xwei,xArgdat(:,:,n)*X(:,n)-xdata(:,n),K);
       xmole = 0;
        for j = 1 : K
            [xPro] = Mltg(xArgdat(:,:,n)*X(:,n)-xdata(:,n),Sig(:,:,j));
           xmole = xmole + xwei(j)*xPro*xArgdat(:,:,n)'*inv(Sig(:,:,j));
        end
       partXn  = xmole/xweisum;
    end