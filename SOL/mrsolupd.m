function [u2, counts] = mrsolupd(u1,rid, nlvl,lvb,bvps,bghs, FI,SI,FE,SE,sz, A1,A2, Mb,pd)
% direct solution of a perturbed problem
% A2*u2 = f - A2*u1 = (A1-A2)*u1;

nbd = length(Mb);
rlvl = floor(log2(rid)) + 1;
nod   = 2^(nlvl-rlvl);
nodtr = 2^(nlvl-1);
r = zeros(size(pd,1),nod);
for il = nod:nod*2-1
    iltr = il + (rid-1)*nod; % index in the global tree
    iatr = iltr - nodtr + 1; % index in the global level
    r(:,il-nod+1) = A1{iatr}*u1(pd(:,iatr))-A2{iatr}*u1(pd(:,iatr));
end

timer = tic;

% refactorization
[FI2,SI2,sz2, counts] = mrfacint(rid,nlvl,lvb,bvps,bghs, A2,Mb);
flp1 = sum(counts.nflops);
stg1 = sum(counts.storage);

% boundary of the node rid
% N: SE{:,:,rid}, M: SI2{:,:,1}
% NM-I;
FIE = cell(nbd,nbd+2);
bds = find(sz2(:,1));
bds = bds(:)';
ne = sum(sz2(bds,1));

for jbd = bds
    for ibd = bds
        if (ibd==jbd)
           FIE{jbd,ibd} = -eye(sz2(jbd,1));
        else
           FIE{jbd,ibd} = zeros(sz2(jbd,1),sz2(ibd,1));
        end
        for kbd = bds
           FIE{jbd,ibd} = FIE{jbd,ibd} + SE{jbd,kbd,rid}*SI2{kbd,ibd,1};
        end
    end
end
flp1 = flp1 + 2*ne^3;
stg1 = stg1 + ne^2;

% block LU
for j = 1:length(bds)
    jbd = bds(j);
    % LU of A11
    [FIE{jbd,nbd+1},FIE{jbd,nbd+2}] = lu(FIE{jbd,jbd});
    FIE{jbd,jbd} = [];
    for i = j+1:length(bds)
        ibd = bds(i);
        % inv(L)*A12
        FIE{jbd,ibd} = FIE{jbd,nbd+1} \ FIE{jbd,ibd};
        % A21*inv(U)
        FIE{ibd,jbd} = FIE{ibd,jbd} / FIE{jbd,nbd+2};
    end
    for i = j+1:length(bds)
        ibd = bds(i);
        for k = j+1:length(bds)
           kbd = bds(k);
           FIE{kbd,ibd} = FIE{kbd,ibd} - FIE{kbd,jbd}*FIE{jbd,ibd};
        end
    end
end
flp1 = flp1 + 2/3*ne^3;


% interior solution
[u2,g,counts] = mrsolint(size(u1),r,rid,nlvl,lvb,bvps,bghs, FIE,SE,FI2,SI2,sz2, Mb,pd);
fprintf('  Fac flops %.2e, storage %.2e; interior solution %.2e\n',flp1,stg1,counts.nflops);
flp1 = flp1 + counts.nflops;
FI2 = [];
SI2 = [];
FIE = [];
t1 = toc(timer);


flp2 = 0.0;
% exterior solution
tic;
mid  = rid;
mlvl = rlvl;
while(mlvl>1)
    pid = floor(mid/2);
    sid = pid*4+1 - mid;
    
    % size for the assembly
    c1 = pid*2;
    c2 = pid*2+1;
    szc = zeros(2,nbd);
    for ibd = 1:nbd
        szc(:,ibd) = [sz(ibd,c1),sz(ibd,c2)];
        if (lvb(ibd,c1) == mlvl)
            szc(1,ibd) = 0;
            ec1 = ibd;
        end
        if (lvb(ibd,c2) == mlvl)
            szc(2,ibd) = 0;
            ec2 = ibd;
        end
    end
    
    bds1 = find(szc(1,:));
    bds2 = find(szc(2,:));
    bds1 = bds1(:)';
    bds2 = bds2(:)';
    ne1  = sum(szc(1,bds1));
    ne2  = sum(szc(2,bds2));
    
    bds = find(sz(:,pid));
    bds = bds(:)';
    for ibd = bds
        g{ibd,mlvl-1} = zeros(sz(ibd,pid),1);
    end
    
    if (mid<sid)
        for ibd = bds1
            g{ibd,mlvl-1}(1:szc(1,ibd)) = g{ibd,mlvl};
        end
        
        %/ KN    K\/-h1\   / -K(Nh1+h2)   \
        %\MKN-I MK/\-h2/ = \-MK(Nh1+h2)+h1/
        tmp = cell(nbd,2);
        for ibd = bds2
            %h1
            tmp{ibd,2} = SI{ibd,ec2,c2}*g{ec1,mlvl};
        end
        for ibd = bds2
            %-h2
            tmp{ibd,1} = zeros(szc(2,ibd),1);
            for kbd = bds1
                tmp{ibd,1} = tmp{ibd,1} - SE{ibd,kbd,pid}(szc(1,ibd)+1:end,1:szc(1,kbd))*g{kbd,mlvl};
            end
            %-h2-Nh1
            for kbd = bds2
                tmp{ibd,1} = tmp{ibd,1} - SE{ibd,kbd,pid}(szc(1,ibd)+1:end,szc(1,kbd)+1:end)*tmp{kbd,2};
            end
        end
        flp2 = flp2 + 2*ne2*(sz(ec2,c2) + ne1 + ne2);
        
        %L\
        for i = 1:length(bds2)
           ibd = bds2(i);
           tmp{ibd,1} = FE{ibd,nbd+1,c2} \ tmp{ibd,1};
           for k = i+1:length(bds2)
               kbd = bds2(k);
               tmp{kbd,1} = tmp{kbd,1} - FE{kbd,ibd,c2} * tmp{ibd,1};
           end
        end
        %U\
        for i =length(bds2):-1:1
           ibd = bds2(i);
           tmp{ibd,1} = FE{ibd,nbd+2,c2} \ tmp{ibd,1};
           for k = i-1:-1:1
               kbd = bds2(k);
               tmp{kbd,1} = tmp{kbd,1} - FE{kbd,ibd,c2} * tmp{ibd,1};
           end
        end
        
        g{ec2,mlvl} = g{ec1,mlvl};
        for ibd = bds2
            g{ibd,mlvl} = tmp{ibd,1};
            for kbd = bds2
                tmp{ibd,2} = tmp{ibd,2} + SI{ibd,kbd,c2}*tmp{kbd,1};
            end
            g{ibd,mlvl-1}(szc(1,ibd)+1:end) = tmp{ibd,2};
        end
        flp2 = flp2 + 4*ne2*ne2;
    else
        for ibd = bds2
            g{ibd,mlvl-1}(szc(1,ibd)+1:end) = g{ibd,mlvl};
        end
        
        %/ KN    K\/-h1\   / -K(Nh1+h2)   \
        %\MKN-I MK/\-h2/ = \-MK(Nh1+h2)+h1/
        
        tmp = cell(nbd,2);
        for ibd = bds1
            %h1
            tmp{ibd,2} = SI{ibd,ec1,c1}*g{ec2,mlvl};
        end
        for ibd = bds1
            %-h2
            tmp{ibd,1} = zeros(szc(1,ibd),1);
            for kbd = bds2
                tmp{ibd,1} = tmp{ibd,1} - SE{ibd,kbd,pid}(1:szc(1,ibd),szc(1,kbd)+1:end)*g{kbd,mlvl};
            end
            %-h2-Nh1
            for kbd = bds1
                tmp{ibd,1} = tmp{ibd,1} - SE{ibd,kbd,pid}(1:szc(1,ibd),1:szc(1,kbd))*tmp{kbd,2};
            end
        end
        flp2 = flp2 + 2*ne1*(sz(ec1,c1) + ne1 + ne2);
        
        %L\
        for i = 1:length(bds1)
           ibd = bds1(i);
           tmp{ibd,1} = FE{ibd,nbd+1,c1} \ tmp{ibd,1};
           for k = i+1:length(bds1)
               kbd = bds1(k);
               tmp{kbd,1} = tmp{kbd,1} - FE{kbd,ibd,c1} * tmp{ibd,1};
           end
        end
        %U\
        for i =length(bds1):-1:1
           ibd = bds1(i);
           tmp{ibd,1} = FE{ibd,nbd+2,c1} \ tmp{ibd,1};
           for k = i-1:-1:1
               kbd = bds1(k);
               tmp{kbd,1} = tmp{kbd,1} - FE{kbd,ibd,c1} * tmp{ibd,1};
           end
        end
        
        g{ec1,mlvl} = g{ec2,mlvl};
        for ibd = bds1
            g{ibd,mlvl} = tmp{ibd,1};
            for kbd = bds1
                tmp{ibd,2} = tmp{ibd,2} + SI{ibd,kbd,c1}*tmp{kbd,1};
            end
            g{ibd,mlvl-1}(1:szc(1,ibd)) = tmp{ibd,2};
        end
        flp2 = flp2 + 4*ne1*ne1;
    end
    [u2,flp] = mrsolbvp(u2,g,sid,nlvl,lvb,bvps,bghs,FI,SI,sz, Mb,pd);
    flp2 = flp2 + flp;
    mid  = pid;
    mlvl = mlvl-1;
end
t2 = toc;

counts = struct('time', [t1,t2],'nflops',[flp1,flp2],'storage',stg1);
end

function [u, flp] = mrsolbvp(u,g,rid,nlvl,lvb,bvps,bghs,FI,SI,sz, Mb,pd)
nbd = 4;
flp = 0.0;
rlvl = floor(log2(rid)) + 1;
nod = 2^(nlvl-rlvl+1)-1; % number of nodes in the subtree


gi = cell(nbd,nod);
for ibd = 1:4
    gi{ibd,1} = g{ibd,rlvl};
end

for lvl = rlvl:nlvl-1
    nod   = 2^(lvl-rlvl);
    for il = nod:nod*2-1
        iltr = il + (rid-1)*nod; % index in the global tree
        c1   = il*2;
        c2   = il*2+1;
        c1tr = iltr*2;
        c2tr = iltr*2+1;
        
        szc = zeros(2,nbd);
        for ibd = 1:nbd
            szc(:,ibd) = [sz(ibd,c1tr),sz(ibd,c2tr)];
            if (lvb(ibd,c1tr) == lvl+1)
                szc(1,ibd) = 0;
                ec1 = ibd;
            end
            if (lvb(ibd,c2tr) == lvl+1)
                szc(2,ibd) = 0;
                ec2 = ibd;
            end
        end
        
        % interface
        ne = sz(ec1,c1tr);
        tmp1 = zeros(ne,1);
        tmp2 = zeros(ne,1);
        bds = find(sz(:,iltr));
        bds = bds(:)';
        for ibd = bds
            if ( szc(1,ibd) > 0 )
                gi{ibd,c1} = gi{ibd,il}(1:szc(1,ibd));
                tmp1 = tmp1 + SI{ec1,ibd,c1tr}*gi{ibd,c1};
                flp = flp + 2*ne*sz(ibd,c1tr);
            end
            if ( szc(2,ibd) > 0 )
                gi{ibd,c2} = gi{ibd,il}(szc(1,ibd)+1:end);
                tmp2 = tmp2 + SI{ec2,ibd,c2tr}*gi{ibd,c2};
                flp = flp + 2*ne*sz(ibd,c2);
            end
            gi{ibd,il} = [];
        end
        %/ KN    K\/-h1\   / -K(Nh1+h2)   \
        %\MKN-I MK/\-h2/ = \-MK(Nh1+h2)+h1/
        tmp2 = FI{2,iltr}\(FI{1,iltr}\(-tmp2-SI{ec2,ec2,c2tr}*tmp1 ));
        gi{ec1,c1} = tmp2;
        gi{ec2,c2} = tmp1 + SI{ec1,ec1,c1tr}*tmp2;
        flp = flp + 2*3*ne^2;
    end
end
nod   = 2^(nlvl-rlvl);
nodtr = 2^(nlvl-1);
for il = nod:nod*2-1
    iltr = il + (rid-1)*nod; % index in the global tree
    iatr = iltr - nodtr + 1; % index in the global level
    bvpi = [2, 2, 2, 2]; % Robin
    for ibd = 1:nbd
        if (lvb(ibd,iltr) == 1)
            % global boundary
            bvpi(ibd) = bvps(ibd);
        end
    end
    tmp = zeros(size(pd,1),1);
    for ibd = 1:nbd
        if (sz(ibd,iltr) == 0); continue; end
        gn = bghs(1,2,bvpi(ibd));
        tmp = tmp + Mb{ibd}*gi{ibd,il}/gn;
    end
    u( pd(FI{4,iltr},iatr) ) = u( pd(FI{4,iltr},iatr) ) + FI{2,iltr}\(FI{1,iltr}\tmp(FI{3,iltr}));
    flp = flp + 2*(nnz(FI{1,iltr}) + nnz(FI{2,iltr}));
end
end
