function[weisum] = weightsumgaupdf( Sig,weig,x,K)%
       % The sum of Gaussian components with weights
        Sw =0;
        for j = 1 :K
            [Pro] = Mltg(x,Sig(:,:,j));
            Sw = Sw + Pro*weig(j);
        end
       weisum = Sw; 
    end