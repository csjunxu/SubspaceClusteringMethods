function [Z,p,q] = NVR3(X,r,lambda,gamma1,eta1,gamma2,eta2,maxiter,toler)

[~,~,n] = size(X);  In = eye(n);    Z = In;

XK = getXK(X);  G = getG(XK);   XK2 = getXK2(X);    G2 = getG(XK2);
funval0 = inf;

for t = 1:maxiter
    
    [H1,H2,H3] = getHs(XK,Z);
    H4 = getH4(XK,Z);
    
    [F1,F2,F3] = getHs(XK2,Z);
    F4 = getH4(XK2,Z);
    
    H = H1 - 2*H2 + H3 + gamma1*G + eta1/2*H4;
    H = (H+H')/2;
    
    F = F1 - 2*F2 + F3 + gamma2*G2 + eta2/2*F4;
    F = (F+F')/2;
    
    p = getp(H,r);    q = getp(F,r);
    
    Xp = circprod_p(X,p);    Xq = circprod_q(X,q);
    
    Wp = Xp'*Xp;     Dp = sum(Wp,2);   Lp = diag(Dp) - Wp;
    Wq = Xq'*Xq;     Dq = sum(Wq,2);   Lq = diag(Dq) - Wq;
    
    Z = lyap(Wp+Wq+lambda/2*In,eta1*Lp+lambda/2*In+eta2*Lq,-Wp-Wq);
    %     Z = sylvester(W+lambda/2*In,eta*L+lambda/2*In,W);
    
    funval = sum(sum((Xp-Xp*Z).^2)) + sum(sum((Xq-Xq*Z).^2)) + lambda*sum(sum(Z.^2)) + gamma1*trace(p'*G*p) + gamma2*trace(q'*G2*q) + eta1*trace(Z*Wp*Z') + eta2*trace(Z*Wq*Z');
    err = abs(funval - funval0);
    
    funval0 = funval;

    if err <= toler;    break;  end
    
end


end

function Y = circprod_p(X,p)
[a,~,n] = size(X);
[~,r] = size(p);
Y = zeros(a*r,n);
for t = 1:n
    y = X(:,:,t)*p;
    Y(:,t) = reshape(y,a*r,1);
end

end

function Y = circprod_q(X,q)
[~,b,n] = size(X);
[~,r] = size(q);
Y = zeros(b*r,n);
for t = 1:n
    y = X(:,:,t)'*q;
    Y(:,t) = reshape(y,b*r,1);
end

end

function H4 = getH4(XK,Z)
distZ = dist(Z).^2;
n = size(XK,4);
H4 = 0;
for t = 1:n
    for s = 1:n
        H4 = H4 + distZ(t,s)*XK(:,:,t,s);
    end
end
end

function G = getG(XK)
G = get_kernelG(XK);

if rank(G) < size(G,1)
    G = G + eye(size(G,1))*(1e-8);
end

G = G^(-1);
end

function [ p ] = getp(H,r)
[up,vp] = eig(H);
[~,ivp] = sort(diag(vp),'ascend');
ivp = ivp(1:r);
yvp = up(:,ivp);
p = yvp;
end


function [H1,H2,H3] = getHs(XK,Z)
[n] = size(XK,3);
H1 = 0; H2 = 0; H3 = 0;
for t = 1:n
    zt = Z(t,:);
    for s = 1:n
        zs = Z(s,:);
        zst = zs*zt';
        
        xst = XK(:,:,s,t);
        
        H2 = H2 + Z(t,s)*xst;
        H3 = H3 + zst*xst;
        if t==s
            H1 = H1 + xst;
        end
    end
end
end

function [XK] = getXK(X)
[~,b,n] = size(X);
XK = zeros(b,b,n,n);
for t = 1:n
    xt = X(:,:,t);
    for s = 1:t
        xs = X(:,:,s);
        xst = xt'*xs;
        XK(:,:,t,s) = xst;
    end
end

for t = 1:n
    for s = t:n
        xst = XK(:,:,s,t);
        XK(:,:,t,s) = xst;
    end
end

end

function [XK] = getXK2(X)
[a,~,n] = size(X);
XK = zeros(a,a,n,n);
for t = 1:n
    xt = X(:,:,t);
    for s = 1:t
        xs = X(:,:,s);
        xst = xt*xs';
        XK(:,:,t,s) = xst;
    end
end

for t = 1:n
    for s = t:n
        xst = XK(:,:,s,t);
        XK(:,:,t,s) = xst;
    end
end

end

function G = get_kernelG(XK)
[n] = size(XK,4);
x = 0;
for t = 1:n
    x = x + XK(:,:,t,t);
end
y = sum(sum(XK,4),3);

G = x/n - y/n/n;

end