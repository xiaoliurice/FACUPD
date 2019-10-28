function [y, flp] = ux(nlvl,x,V,rk, sh,ih)
if(nlvl==1)
    y = V{1,sh,ih}(:,1:rk(sh,ih))*x;
    flp = 2*length(y)*length(x);
    return;
end

flp = 0.0;
vxt = cell(2^(nlvl-1),1);

vxt{1} = x;
nod = 1;
nv = ih;
for lvl = 1:nlvl
    for it=nod:-1:1
        iv=it+nv-1;
        
        vxt{it} = V{1,sh,iv}(:,1:rk(sh,iv))*vxt{it};
        flp = flp + 2*rk(sh,iv)*length(vxt{it});
        
        if(lvl==nlvl); continue; end
        c1t = it*2-1;
        c2t = it*2;
        c1v = iv*2;
        cut = rk(sh,c1v);
        vxt{c2t} = vxt{it}(cut+1:end);
        vxt{c1t} = vxt{it}(1:cut);
    end
    nod = nod*2;
    nv  = nv*2;
end
y=vertcat(vxt{:});
end