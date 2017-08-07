function NcutValue = ComputeNcutValue(W,grp)
K = size(grp,2);
NcutValue = 0;

for ii = 1:K
	idx = grp(:,ii)==1;
	Wg = W(~idx,idx);
	cutA = sum(Wg(:));	
	NcutValue = NcutValue+cutA;
end

end