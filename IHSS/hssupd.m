function [Bo, flp, nnz] = hssupd(nlvl,V,rk,Bi,M,sh,ih)

if (nlvl == 1)
    Bo = Bi + V{1,sh,ih}(:,1:size(M,1)) * M * V{2,sh,ih}(:,1:size(M,2))';
    flp = 2*(size(Bo,1)*numel(M) + numel(Bo)*size(M,2));
    nnz = numel(Bo);
    return;
end

n1d      = 2^nlvl;
Bo{n1d}  = M;
flp      = 0.0;
nnz      = 0.0;
nod = 1;
n1v = ih;
for lvl = 1:nlvl-1
    for ib = nod*2-1:-1:nod
        id  = ib+n1d-nod;
        iv  = ib+n1v-nod;
        T = V{1,sh,iv}(:,1:size(Bo{id},1))*Bo{id}*V{2,sh,iv}(:,1:size(Bo{id},2))';
        flp = flp + 2*(size(T,1)*numel(Bo{id}) + numel(T)*size(Bo{id},2));
        
        c1b = ib*2;
        c2b = ib*2+1;
        c1d = (ib-nod)*2+n1d;
        c2d = c1d+1;
     
        rk1 = rk(sh,iv*2);
        rk2 = rk1;
        
        Bo{c1b} = Bi{c1b} + T(1:rk1,rk2+1:end);
        Bo{c2b} = Bi{c2b} + T(rk1+1:end,1:rk1);
        nnz = nnz + numel(Bo{c1b}) + numel(Bo{c2b});
        Bo{c2d} = T(rk1+1:end,rk2+1:end);
        Bo{c1d} = T(1:rk1,1:rk2);
    end
    nod = nod*2;
    n1v = n1v*2;
end

for id = n1d:n1d+nod-1
    iv = id - n1d + n1v;
    flp = flp + 2*(size(T,1)*numel(Bo{id}) + numel(T)*size(Bo{id},2));
    Bo{id} = Bi{id} + V{1,sh,iv}(:,1:size(Bo{id},1))*Bo{id}*(V{2,sh,iv}(:,1:size(Bo{id},2)))';
    nnz = nnz + numel(Bo{id});
end

end