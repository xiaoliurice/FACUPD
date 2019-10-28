function [u,g,counts] = mrsolhssint(f,r,rid,nlvl,slvl,lvb,bvps,bghs, FIE,SE,V2,rk2,FI,SI,sz,V,rk, lvh,whb,shb,ihb,lhb, M,Mb,pd)
% direct solution of an interior problem (rid root)
% Au = Mf + r;
u = zeros(size(f));
rlvl = floor(log2(rid)) + 1;
flp = 0.0;

nbd = 4;
nod = 2^(nlvl-rlvl+1)-1; % number of nodes in the subtree
gi = cell(nbd,nod);
hi = cell(nbd,nod);

tic;
nod   = 2^(nlvl-rlvl);
nodtr = 2^(nlvl-1);
for il = nod:nod*2-1
    iltr = il + (rid-1)*nod; % index in the global tree
    iatr = iltr - nodtr + 1; % index in the global level
    bvpi = [2, 2, 2, 2];     % Robin
    for ibd = 1:nbd
        if (lvb(ibd,iltr) == 1)
            % global boundary
            bvpi(ibd) = bvps(ibd);
        end
    end
    tmp = M*f(pd(:,iatr)) + r(:,il-nod+1);
    tmp(FI{4,il}) = FI{2,il}\(FI{1,il}\tmp(FI{3,il}));
%    u(pd(:,iatr)) = u(pd(:,iatr)) + tmp;
    flp = flp + 2*(nnz(FI{1,il}) + nnz(FI{2,il}));
    
    for ibd = 1:nbd
        if (sz(ibd,il) == 0); continue; end
        hd = bghs(2,1,bvpi(ibd));
        hn = bghs(2,2,bvpi(ibd));
        gd = bghs(1,1,bvpi(ibd)); 
        gn = bghs(1,2,bvpi(ibd));
        const = hd-gd*hn/gn;
        
        hi{ibd,il} = const * (Mb{ibd}'*tmp);
    end
end


for lvl = nlvl-1:-1:max(rlvl,slvl)
    nod   = 2^(lvl-rlvl);
    for il = nod:nod*2-1
        iltr = il + (rid-1)*nod; % index in the global tree
        
        c1   = il*2;
        c2   = il*2+1;
        c1tr = iltr*2;
        c2tr = iltr*2+1;
        
        szc = zeros(2,nbd);
        for ibd = 1:nbd
            szc(:,ibd) = [sz(ibd,c1),sz(ibd,c2)];
            if (lvb(ibd,c1tr) == lvl+1)
                szc(1,ibd) = 0;
                ec1 = ibd;
            end
            if (lvb(ibd,c2tr) == lvl+1)
                szc(2,ibd) = 0;
                ec2 = ibd;
            end
        end
        %/ KN    K\/-h1\   / -K(Nh1+h2)   \
        %\MKN-I MK/\-h2/ = \-MK(Nh1+h2)+h1/
        % interface
        gi{ec1,c1} = FI{2,il}\(FI{1,il}\ (-hi{ec2,c2}-SI{ec2,ec2,c2}*hi{ec1,c1}) );
        gi{ec2,c2} = hi{ec1,c1} + SI{ec1,ec1,c1}*gi{ec1,c1};

        ne = sz(ec1,c1);
        flp = flp + 2*3*ne^2;

        % g{lvl+1} --> h{lvl}
        bds = find(sz(:,il)); % nonzero blocks
        bds = bds(:)';
        for ibd = bds
            hi{ibd,il} = zeros(sz(ibd,il),1);
            if ( szc(1,ibd) > 0 )
                hi{ibd,il}(1:szc(1,ibd))     = hi{ibd,c1} + SI{ibd,ec1,c1}*gi{ec1,c1};
                flp = flp + 2*sz(ibd,c1)*ne;
            end
            if ( szc(2,ibd) > 0 )
                hi{ibd,il}(szc(1,ibd)+1:end) = hi{ibd,c2} + SI{ibd,ec2,c2}*gi{ec2,c2};
                flp = flp + 2*sz(ibd,c2)*ne;
            end
        end
    end
end

for lvl = slvl-1:-1:rlvl
    nod   = 2^(lvl-rlvl);
    for il = nod:nod*2-1
        iltr = il + (rid-1)*nod; % index in the global tree
        
        c1   = il*2;
        c2   = il*2+1;
        c1tr = iltr*2;
        c2tr = iltr*2+1;
        
        szc = zeros(2,nbd);
        for ibd = 1:nbd
            szc(:,ibd) = [sz(ibd,c1),sz(ibd,c2)];
            if (lvb(ibd,c1tr) == lvl+1)
                szc(1,ibd) = 0;
                ec1 = ibd;
            end
            if (lvb(ibd,c2tr) == lvl+1)
                szc(2,ibd) = 0;
                ec2 = ibd;
            end
        end
        wh = whb(ec1,c1tr);
        lh = lvh(wh);
        [gi{ec1,c1}, gi{ec2,c2}, flp1] = hsslusolint(lh,-hi{ec1,c1},-hi{ec2,c2}, FI{1,il},FI{2,il},FI{3,il},FI{4,il}, V{wh},rk{wh} );
        flp = flp + flp1;
        if(sum(szc(1,:))>0)
            [vg1, flp1] = vtx(lh,gi{ec1,c1}, V{wh},rk{wh},1,1);
            flp = flp + flp1;
        end
        if(sum(szc(2,:))>0)
            [vg2, flp1] = vtx(lh,gi{ec2,c2}, V{wh},rk{wh},2,1);
            flp = flp + flp1;
        end
        
        % g{lvl+1} --> h{lvl}
        bds = find(sz(:,il)); % nonzero blocks
        bds = bds(:)';
        for ibd = bds
            hi{ibd,il} = zeros(sz(ibd,il),1);
            if ( szc(1,ibd) > 0 )
                wh = whb(ibd,c1tr);
                lh = lvh(wh) - lhb(ibd,c1tr) + 1;
                [sg,flp1] = ux(lh,SI{ibd,ec1,c1}*vg1, V{wh},rk{wh},shb(ibd,c1tr),ihb(ibd,c1tr));
                flp = flp + flp1 + 2*numel(SI{ibd,ec1,c1});
                hi{ibd,il}(1:szc(1,ibd))     = hi{ibd,c1} + sg;
            end
            if ( szc(2,ibd) > 0 )
                wh = whb(ibd,c2tr);
                lh = lvh(wh) - lhb(ibd,c2tr) + 1;
                [sg, flp1] = ux(lh,SI{ibd,ec2,c2}*vg2, V{wh},rk{wh},shb(ibd,c2tr),ihb(ibd,c2tr));
                flp = flp + flp1 + 2*numel(SI{ibd,ec1,c1});
                hi{ibd,il}(szc(1,ibd)+1:end) = hi{ibd,c2} + sg;
            end
        end
    end
end

% coupling on the boundary
bds = find(sz(:,1)); % nonzero blocks
bds = bds(:)';
g = cell(nbd,rlvl);
if(rlvl>slvl)
    ne = sum(sz(bds,1));
    for ibd = bds
        gi{ibd,1} = zeros(sz(ibd,1),1);
        for kbd = bds
            gi{ibd,1} =  gi{ibd,1} - SE{ibd,kbd,rid}*hi{kbd,1};
        end
    end
    %L\
    for i = 1:length(bds)
        ibd = bds(i);
        gi{ibd,1} = FIE{ibd,nbd+1} \ gi{ibd,1};
        for k = i+1:length(bds)
            kbd = bds(k);
            gi{kbd,1} = gi{kbd,1} - FIE{kbd,ibd} * gi{ibd,1};
        end
    end
    %U\
    for i = length(bds):-1:1
        ibd = bds(i);
        gi{ibd,1} = FIE{ibd,nbd+2} \ gi{ibd,1};
        for k = i-1:-1:1
            kbd = bds(k);
            gi{kbd,1} = gi{kbd,1} - FIE{kbd,ibd} * gi{ibd,1};
        end
    end
    
    
    for ibd = bds
        g{ibd,rlvl} = hi{ibd,1};
        for kbd = bds
            g{ibd,rlvl} = g{ibd,rlvl} + SI{ibd,kbd,1}*gi{kbd,1};
        end
    end
    flp = flp + 2*3*ne^2;
else
    % put hss solution here
    tb = cell(nbd,1);
    tc = cell(nbd,1);
    for ibd = bds
        wh  = whb(ibd,rid);
        sh  = shb(ibd,rid);
        ih  = ihb(ibd,rid);
        lh  = lvh(wh)-lhb(ibd,rid)+1;
        [tb{ibd}, tc{ibd},flp1] = hsslsolext(lh,-hi{ibd,1},zeros(size(hi{ibd,1})),...
                                  FIE{ibd,nbd+3},V{wh},V2{wh},rk{wh},rk2{wh},ih,sh);
        flp = flp +flp1;
    end
    
    %L\
    ne = 0;
    for i = 1:length(bds)
        ibd = bds(i);
        tc{ibd} = FIE{ibd,nbd+1} \ tc{ibd};
        ne = ne + length(tc{ibd});
        for k = i+1:length(bds)
            kbd = bds(k);
            tc{kbd} = tc{kbd} - FIE{kbd,ibd} * tc{ibd};
        end
    end
    %U\
    for i = length(bds):-1:1
        ibd = bds(i);
        tc{ibd} = FIE{ibd,nbd+2} \ tc{ibd};
        for k = i-1:-1:1
            kbd = bds(k);
            tc{kbd} = tc{kbd} - FIE{kbd,ibd} * tc{ibd};
        end
    end
    flp = flp + 2*ne^2;

    for ibd = bds
        wh  = whb(ibd,rid);
        sh  = shb(ibd,rid);
        ih  = ihb(ibd,rid);
        lh  = lvh(wh)-lhb(ibd,rid)+1;
        [gi{ibd,1}, g{ibd,rlvl},flp1] = hssusolext(lh,tb{ibd},tc{ibd},...
                                        FIE{ibd,nbd+4},V{wh},V2{wh},rk{wh},rk2{wh},ih,sh);
        flp = flp +flp1;
    end
    tb = [];
    tc = [];
end

for lvl = rlvl:slvl-1
    nod   = 2^(lvl-rlvl);
    for il = nod:nod*2-1
        iltr = il + (rid-1)*nod; % index in the global tree
        c1   = il*2;
        c2   = il*2+1;
        c1tr = iltr*2;
        c2tr = iltr*2+1;
        
        szc = zeros(2,nbd);
        for ibd = 1:nbd
            szc(:,ibd) = [sz(ibd,c1),sz(ibd,c2)];
            if (lvb(ibd,c1tr) == lvl+1)
                szc(1,ibd) = 0;
                ec1 = ibd;
            end
            if (lvb(ibd,c2tr) == lvl+1)
                szc(2,ibd) = 0;
                ec2 = ibd;
            end
        end
        
        wh = whb(ec1,c1tr);
        lh = lvh(wh);
        if(sum(szc(1,:)) > 0)
            tmp1 = zeros(rk{wh}(1,1),1);
        end
        if(sum(szc(2,:)) > 0)
            tmp2 = zeros(rk{wh}(2,1),1);
        end
        
        bds = find(sz(:,il)); % nonzero blocks
        bds = bds(:)';
        for ibd = bds
            if(szc(1,ibd) > 0)
                gi{ibd,c1} = gi{ibd,il}(1:szc(1,ibd));
                whc = whb(ibd,c1tr);
                [vg1, flp1] = vtx(lvh(whc)-lhb(ibd,c1tr)+1,gi{ibd,c1}, V{whc},rk{whc},shb(ibd,c1tr),ihb(ibd,c1tr));
                tmp1 = tmp1 + SI{ec1,ibd,c1}*vg1;
                flp = flp + flp1 + 2*length(tmp1)*length(vg1);
            end
            if(szc(2,ibd) > 0)
                gi{ibd,c2} = gi{ibd,il}(szc(1,ibd)+1:end);
                whc = whb(ibd,c2tr);
                [vg2, flp1] = vtx(lvh(whc)-lhb(ibd,c2tr)+1,gi{ibd,c2}, V{whc},rk{whc},shb(ibd,c2tr),ihb(ibd,c2tr));
                tmp2 = tmp2 + SI{ec2,ibd,c2}*vg2;
                flp = flp + flp1 + 2*length(tmp2)*length(vg2);
            end
        end
        
        if(sum(szc(1,:)) > 0)
            [tx1, flp1] = ux(lh,tmp1, V{wh},rk{wh},1,1);
            flp = flp + flp1;
        else
            tx1 = zeros(sz(ec1,c1),1);
        end
        
        if(sum(szc(2,:)) > 0)
            [tx2, flp1] = ux(lh,tmp2, V{wh},rk{wh},2,1);
            flp = flp + flp1;
        else
            tx2 = zeros(sz(ec2,c2),1);
        end

        [tx1,tx2, flp1] = hsslusolint(lh,-tx1,-tx2, FI{1,il},FI{2,il},FI{3,il},FI{4,il}, V{wh},rk{wh} );
        flp = flp + flp1;
        
        gi{ec1,c1} = gi{ec1,c1} + tx1;
        gi{ec2,c2} = gi{ec2,c2} + tx2;
    end
end


for lvl = max(rlvl,slvl):nlvl-1
    nod   = 2^(lvl-rlvl);
    for il = nod:nod*2-1
        iltr = il + (rid-1)*nod; % index in the global tree
        c1   = il*2;
        c2   = il*2+1;
        c1tr = iltr*2;
        c2tr = iltr*2+1;
        
        szc = zeros(2,nbd);
        for ibd = 1:nbd
            szc(:,ibd) = [sz(ibd,c1),sz(ibd,c2)];
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
        ne = sz(ec1,c1);
        hi{ec1,c1} = zeros(ne,1);
        hi{ec2,c2} = zeros(ne,1);
        bds = find(sz(:,il));
        bds = bds(:)';
        
        for ibd = bds          
            if ( szc(1,ibd) > 0 )
                gi{ibd,c1} = gi{ibd,il}(1:szc(1,ibd));
                hi{ec1,c1} = hi{ec1,c1} + SI{ec1,ibd,c1}*gi{ibd,c1};
                flp = flp + 2*ne*sz(ibd,c1);
            end
            if ( szc(2,ibd) > 0 )
                gi{ibd,c2} = gi{ibd,il}(szc(1,ibd)+1:end);
                hi{ec2,c2} = hi{ec2,c2} + SI{ec2,ibd,c2}*gi{ibd,c2};
                flp = flp + 2*ne*sz(ibd,c2);
            end
        end
        %/ KN    K\/-h1\   / -K(Nh1+h2)   \
        %\MKN-I MK/\-h2/ = \-MK(Nh1+h2)+h1/
        hi{ec2,c2} = FI{2,il}\(FI{1,il}\(-hi{ec2,c2}-SI{ec2,ec2,c2}*hi{ec1,c1} ));
        gi{ec1,c1} = gi{ec1,c1} + hi{ec2,c2};
        gi{ec2,c2} = gi{ec2,c2} + hi{ec1,c1} + SI{ec1,ec1,c1}*hi{ec2,c2};
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
    tmp = M*f(pd(:,iatr)) + r(:,il-nod+1);
    for ibd = 1:nbd
        if (sz(ibd,il) == 0); continue; end
        gn = bghs(1,2,bvpi(ibd));
        tmp = tmp + Mb{ibd}*gi{ibd,il}/gn;
    end
    tmp(FI{4,il}) = FI{2,il}\(FI{1,il}\tmp(FI{3,il}));
    flp = flp + 2*(nnz(FI{1,il}) + nnz(FI{2,il}));
    u(pd(:,iatr)) = tmp;
end
t = toc;
counts = struct('time', t,'nflops',flp);
end