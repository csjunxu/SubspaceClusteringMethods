function matrank=modelselection(sv,k)
for(r=1:length(sv)-1)
    allr(r)=sv(r+1)^2/sum(sv(1:r).^2)+k*r;
end
[minr,matrank]=min(allr);
