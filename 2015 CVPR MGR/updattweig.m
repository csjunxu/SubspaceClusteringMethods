% Updating the pth weight
function[newweip] = updattweig(wdatam,owei,p,K,X,Sig)
      [pgaml] = gamij(wdatam,Sig,p,X,K,owei);
      rwdatc = size(wdatam,2);
      newweip = (owei(p)/rwdatc)*pgaml;
      end