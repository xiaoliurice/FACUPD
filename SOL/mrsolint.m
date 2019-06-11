function [u,g,counts] = mrsolint(sizeu,r,rid,nlvl,lvb,bvps,bghs, FIE,SE,FI,SI,sz, Mb,pd)
% direct solution of an interior problem (rid root)
% Au = r;
u = zeros(sizeu);

rlvl = floor(log2(rid)) + 1;
flp = 0.0;

nbd = 4;
nod = 2^(nlvl-rlvl+1)-1; % number of nodes in the subtree
gi = cell(nbd,nod);
hi = cell(nbd,nod);

tic;
nod   = 2^(nlvl-rlvl);
nodtr = 2^(nlvl-1);
tmp = zeros(size(pd,1),1);
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
    tmp(FI{4,il}) = FI{2,il}\(FI{1,il}\r(FI{3,il},il-nod+1));
    u(pd(:,iatr)) = u(pd(:,iatr)) + tmp;
    flp = flp + 2*(nnz(FI{1,il}) + nnz(FI{2,il}));
    
    for ibd = 1:nbd
        if (sz(ibd,il) == 0);  continue; end
        hd = bghs(2,1,bvpi(ibd));
        hn = bghs(2,2,bvpi(ibd));
        gd = bghs(1,1,bvpi(ibd)); 
        gn = bghs(1,2,bvpi(ibd));
        const = hd-gd*hn/gn;
        
        hi{ibd,il} = const * (Mb{ibd}'*tmp);
    end
end


for lvl = nlvl-1:-1:rlvl
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
            hi{ibd,c1} = [];
            hi{ibd,c2} = [];
        end
    end
end


% coupling on the boundary
bds = find(sz(:,1)); % nonzero blocks
bds = bds(:)';
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

g = cell(nbd,rlvl);
for ibd = bds
    g{ibd,rlvl} = hi{ibd,1};
    for kbd = bds
        g{ibd,rlvl} = g{ibd,rlvl} + SI{ibd,kbd,1}*gi{kbd,1};
    end
end
flp = flp + 2*3*ne^2;

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
        tmp1 = zeros(ne,1);
        tmp2 = zeros(ne,1);
        bds = find(sz(:,il));
        bds = bds(:)';
        
        for ibd = bds
            if ( szc(1,ibd) > 0 )
                gi{ibd,c1} = gi{ibd,il}(1:szc(1,ibd));
                tmp1 = tmp1 + SI{ec1,ibd,c1}*gi{ibd,c1};
                flp = flp + 2*ne*sz(ibd,c1);
            end
            if ( szc(2,ibd) > 0 )
                gi{ibd,c2} = gi{ibd,il}(szc(1,ibd)+1:end);
                tmp2 = tmp2 + SI{ec2,ibd,c2}*gi{ibd,c2};
                flp = flp + 2*ne*sz(ibd,c2);
            end
            gi{ibd,il} = [];
        end        %for ee = ei:ei+bd(ibd)-1
        %/ KN    K\/-h1\   / -K(Nh1+h2)   \
        %\MKN-I MK/\-h2/ = \-MK(Nh1+h2)+h1/
        tmp2 = FI{2,il}\(FI{1,il}\(-tmp2-SI{ec2,ec2,c2}*tmp1 ));
        gi{ec1,c1} = gi{ec1,c1} + tmp2;
        gi{ec2,c2} = gi{ec2,c2} + tmp1 + SI{ec1,ec1,c1}*tmp2;
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
        if (sz(ibd,il) == 0); continue; end
        gn = bghs(1,2,bvpi(ibd));
        tmp = tmp + Mb{ibd}*gi{ibd,il}/gn;
    end
    u(pd(FI{4,il},iatr)) = u(pd(FI{4,il},iatr)) + FI{2,il}\(FI{1,il}\tmp(FI{3,il}));
    flp = flp + 2*(nnz(FI{1,il}) + nnz(FI{2,il}));
end
t = toc;
counts = struct('time', t,'nflops',flp);
end