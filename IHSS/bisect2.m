function [nsub, lvb, clf, lvh, whb,shb,ihb,lhb,nre] = bisect2(B0, nlvl, slvl)
% nsub:  number of blocks on each direction
% lvb:   level of boundary
% clf:   coordinate of leaf

% HSS related
% lvh:   level of H (differnt length)
% whb:   which    H for boundary
% shb:   side in  H for boundary
% ihb:   index in H for boundary
% lhb:   level of H for boundary

nsub   = ones(size(B0));

clf = ones(2,2^(nlvl-1));
lvb = ones(4,2^nlvl-1);
whb = zeros(4,max(2^slvl-1,1));
shb = whb;
ihb = whb;
lhb = whb;

nn = 1;
oo = 0;
nh = 0;
for lvl = 2: nlvl
    [~,id] = max(B0./nsub);
    nsub(id) = nsub(id)*2;
    no = nn+oo;
    clf(: ,1 : nn*2  ) = kron(clf(:,1:nn),[1,1]);
    clf(id,1:2:nn*2-1) = clf(id,1:2:nn*2-1)*2-1;
    clf(id,2:2:nn*2  ) = clf(id,1:2:nn*2-1)+1;

    lvb(:     ,no+1  :no+nn*2  ) = kron(lvb(:,oo+1:oo+nn),[1,1]);  
    lvb(id*2  ,no+1:2:no+nn*2-1) = lvl;
    lvb(id*2-1,no+2:2:no+nn*2  ) = lvl;
    if (lvl<=slvl)

        whb(:     ,no+1  :no+nn*2  ) = kron(whb(:,oo+1:oo+nn),[1,1]);
        whb(id*2  ,no+1:2:no+nn*2-1) = nh+1:nh+nn;
        whb(id*2-1,no+2:2:no+nn*2  ) = nh+1:nh+nn;
        
        shb(:     ,no+1  :no+nn*2  ) = kron(shb(:,oo+1:oo+nn),[1,1]);
        shb(id*2  ,no+1:2:no+nn*2-1) = 1;
        shb(id*2-1,no+2:2:no+nn*2  ) = 2;
        
        lhb(:     ,no+1  :no+nn*2  ) = kron(lhb(:,oo+1:oo+nn),[1,1])+1;  
        lhb(id*2-1:id*2,no+1:no+nn*2)= lhb(id*2-1:id*2,no+1:no+nn*2)-1;
        lhb(id*2  ,no+1:2:no+nn*2-1) = 1;
        lhb(id*2-1,no+2:2:no+nn*2  ) = 1;
        
        ihb(:     ,no+1  :no+nn*2  ) = kron(ihb(:,oo+1:oo+nn),[1,1]);
        ihb(id*2  ,no+1:2:no+nn*2-1) = 1;
        ihb(id*2-1,no+2:2:no+nn*2  ) = 1;
        for dm = [1:id*2-2, id*2+1:length(B0)*2]
            ihb(dm,no+1:2:no+nn*2-1) = ihb(dm,no+1:2:no+nn*2-1).*2;
            ihb(dm,no+2:2:no+nn*2  ) = ihb(dm,no+2:2:no+nn*2  ).*2+sign( ihb(dm,no+2:2:no+nn*2) );
        end
        
        nh = nh + nn;
    end
    
    oo = no;
    nn = nn*2;
end
nre = [0,0];
lvh = ones(1,nh);
for j = 1: numel(lhb)
    k = whb(j);
    if (k>0)
       lvh(k) = max(lhb(j), lvh(k));
    end
end
lhb = lhb.*sign(ihb);

for j = 1: numel(whb)
    k = whb(j);
    if(k>0 && lvb(j) <= 3 )
        nre(1) = nre(1) + 2^(lvh(k)-lhb(j)+1)-1;
    end
end
nre(2) = sum(2.^lvh(1:min(3,end))-1)*2;

end
