% find the qth covariance matrix
function[Sigmq] = findcovj(fcodat,Sig,K,X,Weicov,q)
        [gamq] = gamij(fcodat,Sig,q,X,K,Weicov);
        covc = size(fcodat,2);
        Sigmq = 0;
         [Arfcodat] = avoidze(fcodat);
        for j = 1 : covc
            [Prosq] = Mltg(Arfcodat(:,:,j)*X(:,j)-fcodat(:,j),Sig(:,:,q));
            [weisumq] = weightsumgaupdf(Sig,Weicov,Arfcodat(:,:,j)*X(:,j)-fcodat(:,j),K);
            Sigfenzi = Prosq*(Arfcodat(:,:,j)*X(:,j)-fcodat(:,j))*(Arfcodat(:,:,j)*X(:,j)-fcodat(:,j))';
            Sigmq = Sigmq + Sigfenzi/weisumq;
        end
        Sigmq = (1/gamq)*Sigmq + 1E-5.*eye(size(Sigmq)) ;
    end% 