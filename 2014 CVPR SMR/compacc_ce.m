function [err,idx_wrong] = compacc_ce(idx,gnd)

% This function computes clustering error (CE)

gnd = compress_labels(gnd);
idx = compress_labels(idx);

[res_lgc, temp] = bestMap(gnd,idx);
if size(res_lgc,1)~=size(gnd,1)
    gnd = gnd';
end
idx_wrong = find(gnd ~= res_lgc);
acc = 1-length(idx_wrong)/length(gnd);

err = 1-acc;


function [gnd_dst,igids] = compress_labels(gnd)

gids = unique(gnd)+1;
igids = zeros(1,max(gids));
for i=1:length(gids)
    igids(gids(i)) = i;
end

gnd_dst = igids(gnd+1);
