function [x1,x2,flp] = hsslusolint(nlvl,b1,b2,L,U,L2,U2,V,rk)
if(nlvl==1)
    x1 =  U\ (L\ (U2*b1 + b2) );
    x2 = -b1 + L2*x1;
    flp = length(b1)^2 * 6;
    return;
end
flp = 0.0;

nod = 2^(nlvl-1);
n1d = 2^nlvl;
ob1 = 0;
ob2 = 0;
tb  = cell(2,n1d-1);
for lvl=nlvl:-1:2
    for iv = nod:nod*2-1
        if(lvl==nlvl)
            szt = [size(V{1,1,iv},1),size(V{1,2,iv},1)];
            tb{1,iv} = V{1,1,iv}'*b1(ob1+1:ob1+szt(1));
            ob1 = ob1 + szt(1);
            tb{2,iv} = V{1,2,iv}'*b2(ob2+1:ob2+szt(2));
            ob2 = ob2 + szt(2);
            flp = flp + 4 * szt(1)^2;
        else
            c1v = iv*2;
            c2v = iv*2+1;
            
            tb{1,iv} = V{1,1,iv}'*[tb{1,c1v}(1:rk(1,c1v));tb{1,c2v}(1:rk(1,c2v))];
            tb{2,iv} = V{1,2,iv}'*[tb{2,c1v}(1:rk(2,c1v));tb{2,c2v}(1:rk(2,c2v))];
            szt = [length(tb{1,iv}),length(tb{2,iv})];
            flp = flp + 2*(szt(1)^2+szt(2)^2);
        end
        
        rkt = rk(:,iv);
        tx1 = L{iv} \ [tb{1,iv}(rkt(1)+1:end);tb{2,iv}(rkt(2)+1:end)];
        tx2 = [tb{1,iv}(1:rkt(1));tb{2,iv}(1:rkt(2))] - L2{iv}*tx1;
        flp = flp + length(tx1)*(length(tx1) + 2*length(tx2));
        
        tb{1,iv} = [tx2(1:rkt(1));     tx1(1:szt(1)-rkt(1))];
        tb{2,iv} = [tx2(rkt(1)+1:end); tx1(szt(1)-rkt(1)+1:end)];
    end
    nod = nod/2;
end

iv = 1;
c1v = iv*2;
c2v = iv*2+1;

tb{1,iv} = [tb{1,c1v}(1:rk(1,c1v));tb{1,c2v}(1:rk(1,c2v))];
tb{2,iv} = [tb{2,c1v}(1:rk(2,c1v));tb{2,c2v}(1:rk(2,c2v))];
cut = length(tb{1,iv});
tx1 = U{iv}\ (L{iv}\ [tb{1,iv};tb{2,iv}]);
flp = flp + 2*length(tx1)^2;
tb{1,c1v}(1:rk(1,c1v)) = tx1(1:rk(1,c1v));
tb{1,c2v}(1:rk(1,c2v)) = tx1(rk(1,c1v)+1:cut);
tb{2,c1v}(1:rk(2,c1v)) = tx1(cut+1:cut+rk(2,c1v));
tb{2,c2v}(1:rk(2,c2v)) = tx1(cut+rk(2,c1v)+1:end);

x1 = zeros(size(b1));
x2 = zeros(size(b2));
ob1 = 0;
ob2 = 0;
nod = 2;
for lvl=2:nlvl
    for iv = nod:nod*2-1
        rkt = rk(:,iv);
        szt = [length(tb{1,iv}),length(tb{2,iv})];
        tx1 = U{iv} \([tb{1,iv}(rkt(1)+1:end);    tb{2,iv}(rkt(2)+1:end)] ...
             -U2{iv}* [tb{1,iv}(1:rkt(1));tb{2,iv}(1:rkt(2))] );
        flp = flp + length(tx1)*( length(tx1) + 2*(rkt(1)+rkt(2)) );
        tb{1,iv}(rkt(1)+1:end) = tx1(1:szt(1)-rkt(1));
        tb{2,iv}(rkt(2)+1:end) = tx1(szt(1)-rkt(1)+1:end);
        
        if(lvl==nlvl)
            szt = [size(V{2,1,iv},1),size(V{2,2,iv},1)];
            x1(ob1+1:ob1+szt(1)) = V{2,1,iv}*tb{1,iv};
            ob1 = ob1 + szt(1);
            x2(ob2+1:ob2+szt(2)) = V{2,2,iv}*tb{2,iv};
            ob2 = ob2 + szt(2);
            flp = flp + 4 * szt(1)^2;
        else
            c1v = iv*2;
            c2v = iv*2+1;
            
            szt = [length(tb{1,iv}),length(tb{2,iv})];
            tb{1,iv} = V{2,1,iv}*tb{1,iv};
            tb{2,iv} = V{2,2,iv}*tb{2,iv};
            flp = flp + 2*(szt(1)^2+szt(2)^2);
            
            tb{1,c1v}(1:rk(1,c1v)) = tb{1,iv}(1:rk(1,c1v));
            tb{1,c2v}(1:rk(1,c2v)) = tb{1,iv}(rk(1,c1v)+1:end);
            tb{2,c1v}(1:rk(2,c1v)) = tb{2,iv}(1:rk(2,c1v));
            tb{2,c2v}(1:rk(2,c2v)) = tb{2,iv}(rk(2,c1v)+1:end);
            
        end
    end
    nod = nod*2;
end
end