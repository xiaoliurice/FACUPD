function [y, flp] = vtx(nlvl,x,V,rk, sh,ih)
if(nlvl==1)
    y = V{2,sh,ih}(:,1:rk(sh,ih))'*x;
    flp = 2*length(y)*length(x);
    return;
end

flp = 0.0;
nod = 2^(nlvl-1);
nv  = ih*nod;
vxt = cell(nod,1);

ox = 0;
for lvl = nlvl:-1:1
    for it = 1:nod
        iv = it-1+nv;
        if(lvl == nlvl)
            szt = size(V{2,sh,iv},1);
            vxt{it} = x(ox+1:ox+szt);
            ox = ox + szt;
        else
            c1t = it*2-1;
            c2t = it*2;
            vxt{it} = [vxt{c1t};vxt{c2t}];
        end
        vxt{it} = V{2,sh,iv}(:,1:rk(sh,iv))'*vxt{it};
        flp = flp + 2*rk(sh,iv)*length(vxt{it});
    end
    nod = nod/2;
    nv  = nv/2;
end
y = vxt{1};
end