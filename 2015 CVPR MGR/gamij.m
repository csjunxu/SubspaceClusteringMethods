 function  [gaml] = gamij(gdat,Sig,l,X,K,Weig)
        gaml = 0;
        [rgdat, cgdat] = size(gdat);
        [Argdat] = avoidze(gdat);
        for j = 1 : cgdat
            [Prog] = Mltg(Argdat(:,:,j)*X(:,j)-gdat(:,j),Sig(:,:,l));
            [weisumg] = weightsumgaupdf(Sig,Weig,Argdat(:,:,j)*X(:,j)-gdat(:,j),K);
            gaml = gaml + Prog/weisumg;
        end
  end