function [FI,SI,sz, V,rk, counts,c1,c2,ec1,ec2] = mrfachssint(rid, nlvl,slvl,tol,lvh,whb,shb,ihb,lhb, lvb,bvps,bghs, A,Mb)
% construct interior factors
% SI: Schur complement as interior boundary maps
% FI: LU factors
% rid: ID of root  (1,1 for root)

rlvl = floor(log2(rid)) + 1;

nbd = 4;
nod = 2^(nlvl-rlvl+1)-1; % number of nodes in the subtree
sz = zeros(nbd,nod);  % size of each BD components

% interior
FI = cell(4,nod); % L,U,p,q;
SI = cell(nbd,nbd,nod);

if (slvl>=rlvl)
    nh = size(lvh,2);
    V  = cell(1,nh);
    rk = cell(1,nh);
    for il = 1:nh
        V{il}  =  cell(2,2, 2^lvh(il)-1);
        rk{il} = zeros(2, 2^lvh(il)-1);
    end
else
    V=0; rk=0;
end

flp1 = 0.0;
flp2 = 0.0;
stg1 = 0.0;
stg2 = 0.0;

%%%%%%%%%%%%%% boundary maps from FEM %%%%%%%%%%%%%%%%%%%%%%%%%
tic;
nod   = 2^(nlvl-rlvl);
nodtr = 2^(nlvl-1);
% local index
for il = nod:nod*2-1
    iltr = il + (rid-1)*nod; % index in the global tree
    iatr = iltr - nodtr + 1; % index in the global level
    bvpi = [2, 2, 2, 2]; % Robin
    
    for ibd = 1:nbd
        if (lvb(ibd,iltr) == 1)
            % global boundary
            bvpi(ibd) = bvps(ibd);
        else
            sz(ibd,il) = size(Mb{ibd},2);
        end
    end
    
    % factorization L U P Q
    [FI{1,il}, FI{2,il}, FI{3,il}, FI{4,il}] = lu(A{iatr},'vector');
    flp1 = flp1 + (2+2*sqrt(2))*(2/3)*(sqrt(nnz(A{iatr})))^3;
    stg1 = stg1 + nnz(FI{1,il}) + nnz(FI{2,il});
    
    for jbd = 1:nbd
        if (sz(jbd,il) == 0); continue; end % skip this portion
        
        gnj = bghs(1,2,bvpi(jbd));
        hnj = bghs(2,2,bvpi(jbd));
        %A\M{2}*I
        tmp = full(Mb{jbd}(FI{3,il},:)/gnj);
        tmp(FI{4,il},:) = FI{2,il}\(FI{1,il}\tmp);
        
        for ibd = 1:nbd
            if (sz(ibd,il) == 0); continue; end
            gd = bghs(1,1,bvpi(ibd));
            gn = bghs(1,2,bvpi(ibd));
            hd = bghs(2,1,bvpi(ibd));
            hn = bghs(2,2,bvpi(ibd));
            const = hd-gd*hn/gn;
            SI{ibd,jbd,il} = const*(Mb{ibd}'*tmp);
        end
        
        SI{jbd,jbd,il} = SI{jbd,jbd,il} + eye(sz(jbd,il),sz(jbd,il))*(hnj/gnj);
        flp1 = flp1 + 2*(nnz(FI{1,il}) + nnz(FI{2,il}))*sz(jbd,il);
    end
    
end
t1=toc;

%%%%%%%%%%%%% merge uncompressed levels %%%%%%%%%%%%%%%%%
tic;
for lvl = nlvl-1:-1:max(slvl,rlvl)
    nod   = 2^(lvl-rlvl);
    for il = nod:nod*2-1
        iltr = il + (rid-1)*nod; % index in the global tree
        
        c1   = il*2;
        c2   = il*2+1;
        c1tr = iltr*2;
        c2tr = iltr*2+1;
        
        % size for the assembly
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
        sz(:,il) = sum(szc);
        
        %  / M,-I\-1
        %  \-I, N/
        %  /  I\/-I         \/I  \
        % =\I M/\  inv(NM-I)/\N I/
        %  /  I\/-I  \ / KN    K\
        % =\I M/\KN K/=\MKN-I MK/
        ne = sz(ec1,c1);
        [FI{1,il},FI{2,il}] = lu(SI{ec2,ec2,c2}*SI{ec1,ec1,c1} - eye(ne));
        stg2 = stg2 + ne^2*3;
        flp2 = flp2 + (2/3+2)*ne^3;
        
        bds = find(sz(:,il)); % nonzero blocks
        bds = bds(:)';
        for jbd = bds
            tmp = cell(2,2);
            stg2 = stg2 + ne*(szc(1,jbd)+szc(2,jbd))*2; % SIN, SNI
            if ( szc(1,jbd) > 0 )
                % inv2x2 * [I;0]
                tmp{1,1} = FI{2,il} \ (FI{1,il} \ (SI{ec2,ec2,c2}*SI{ec1,jbd,c1}));
                tmp{2,1} = SI{ec1,ec1,c1}*tmp{1,1}-SI{ec1,jbd,c1};
                flp2 = flp2 + szc(1,jbd)*(6*ne^2);
            end
            if ( szc(2,jbd) > 0 )
                % inv2x2 * [0;I]
                tmp{1,2} = FI{2,il} \ (FI{1,il} \ SI{ec2,jbd,c2});
                tmp{2,2} = SI{ec1,ec1,c1}*tmp{1,2};
                flp2 = flp2 + szc(2,jbd)*(4*ne^2);
            end
                
            for ibd = bds
                SI{ibd,jbd,il} = zeros(sz(ibd,il),sz(jbd,il));
                stg2 = stg2 + sz(ibd,il)*sz(jbd,il);
                flp2 = flp2 + 2*ne*sz(ibd,il)*sz(jbd,il);
                
                % c1, c1
                if (szc(1,ibd)*szc(1,jbd) > 0)
                    SI{ibd,jbd,il}(1:szc(1,ibd),1:szc(1,jbd))...
                   =SI{ibd,jbd,c1} - SI{ibd,ec1,c1} * tmp{1,1};
                end
                % c1, c2
                if (szc(1,ibd)*szc(2,jbd) > 0)
                    SI{ibd,jbd,il}(1:szc(1,ibd),szc(1,jbd)+1:end)...
                  =-SI{ibd,ec1,c1} * tmp{1,2};
                end
                % c2, c1
                if (szc(2,ibd)*szc(1,jbd) > 0)
                   SI{ibd,jbd,il}(szc(1,ibd)+1:end,1:szc(1,jbd)) ...
                  =-SI{ibd,ec2,c2} * tmp{2,1};
                end
                
                % c2, c2
                if (szc(2,ibd)*szc(2,jbd) > 0)
                    SI{ibd,jbd,il}(szc(1,ibd)+1:end,szc(1,jbd)+1:end)...
                   =SI{ibd,jbd,c2} - SI{ibd,ec2,c2} * tmp{2,2};
                end
            end
        end
        
    end
end


%%%%%%%%%%%%%%%%% 1st level of low-rank compression %%%%%%%%%%%%%%%%%
if (slvl>=rlvl)
    nod = 2^(slvl-rlvl);
    for il = nod:nod*2-1
        iltr = il + (rid-1)*nod; % index in the global tree
        totsz = sum(sz(:,il));
        for jbd = 1:nbd
            if (sz(jbd,il) * (totsz - sz(jbd,il)) == 0 ); continue; end
            wh = whb(jbd,iltr);
            sh = shb(jbd,iltr);
            ih = ihb(jbd,iltr);
            % column compression
            [U1,sig,V1]=svd(vertcat(SI{[1:jbd-1,jbd+1:nbd],jbd,il}));
            flp2 = flp2 + 2*numel(V1)*size(U1,1);
            sig = diag(sig);
            if (sig(end)<=sig(1)*tol)
                rkt = max(find(sig<=sig(1)*tol,1)-1, 1);
            else
                rkt = length(sig);
            end
            rk{wh}(sh,ih) = rkt;
            V{wh}{2,sh,ih} = V1;
            stg2 = stg2 + sz(jbd,il)*rkt;
            for ibd = [1:jbd-1,jbd+1:nbd]
                if (sz(ibd,il) == 0); continue; end
                SI{ibd,jbd,il} = SI{ibd,jbd,il}*V1(:,1:rkt);
            end
        end
        
        % row compression
        for ibd = 1:nbd
            if ( sz(ibd,il) * (totsz - sz(ibd,il)) == 0 ); continue; end
            wh = whb(ibd,iltr);
            sh = shb(ibd,iltr);
            ih = ihb(ibd,iltr);
            [U1,sig,V1]=svd(horzcat(SI{ibd,[1:ibd-1,ibd+1:nbd],il}));
            flp2 = flp2 + 2*numel(U1)*size(V1,1);
            sig = diag(sig);
           
            rkt = rk{wh}(sh,ih);
            V{wh}{1,sh,ih} = U1;
            stg2 = stg2 + sz(ibd,il)*rkt;
            for jbd = [1:ibd-1, ibd+1:nbd]
                if (sz(jbd,il) == 0); continue; end
                SI{ibd,jbd,il} = U1(:,1:rkt)'*SI{ibd,jbd,il};
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
            sz(:,il) = sum(szc);
           
            wh = whb(ec1,c1tr);
            [FI{1,il},FI{2,il},FI{3,il},FI{4,il}, C, flp,stg] = ...
            hssluint(lvh(wh),V{wh},rk{wh},SI{ec1,ec1,c1},SI{ec2,ec2,c2});
            
            stg2 = stg2 + stg;
            flp2 = flp2 + flp;

            bds = find(sz(:,il)); % nonzero blocks
            bds = bds(:)';
            for jbd = bds
                % diagonal
                wh = whb(jbd,iltr);
                sh = shb(jbd,iltr);
                ih = ihb(jbd,iltr);
                lh = lvh(wh)-lhb(jbd,iltr)+1;
                
                ncol1 = min( rk{wh}(shb(jbd,c1tr),ihb(jbd,c1tr)), szc(1,jbd) );
                ncol2 = min( rk{wh}(shb(jbd,c2tr),ihb(jbd,c2tr)), szc(2,jbd) );
                tmp = cell(2);
                flp = 2*numel(C{1,2})*ncol2;
                if(flp > 0)
                    tmp{1,2} = C{1,2}*SI{ec2,jbd,c2};
                    flp2 = flp2 + flp;
                end
                flp = 2*numel(C{2,1})*ncol1;
                if(flp > 0)
                    tmp{2,1} = C{2,1}*SI{ec1,jbd,c1};
                    flp2 = flp2 + flp;
                end
                
                if (ncol1>0)
                    tmp{1,1} = C{1,1}*SI{ec1,jbd,c1};
                    if(ncol2>0)
                        %c1,c2
                        tmp{2,2} = C{2,2}*SI{ec2,jbd,c2};
                        [SI{jbd,jbd,il}, flp] = ...
                        hssupdint2(lh-1,V{wh},rk{wh},SI{jbd,jbd,c1},SI{jbd,jbd,c2},...
                        {-SI{jbd,ec1,c1}*tmp{1,1},-SI{jbd,ec1,c1}*tmp{1,2};...
                         -SI{jbd,ec2,c2}*tmp{2,1},-SI{jbd,ec2,c2}*tmp{2,2}},sh,ih);
                     
                        flp2 = flp2 + flp ...
                        +2*ncol1*(numel(C{1,1})+numel(tmp{1,1})+numel(tmp{1,2}))...
                        +2*ncol2*(numel(C{2,2})+numel(tmp{2,1})+numel(tmp{2,2}));
                    else
                        %c1
                        [SI{jbd,jbd,il}, flp] = ...
                        hssupd(lh,V{wh},rk{wh},SI{jbd,jbd,c1},-SI{jbd,ec1,c1}*tmp{1,1},sh,ih);
                        
                        flp2 = flp2 + flp + 2*ncol1*(numel(C{1,1})+numel(tmp{1,1}));
                    end
                elseif (ncol2>0)
                    %c2
                    tmp{2,2} = C{2,2}*SI{ec2,jbd,c2};
                    
                    [SI{jbd,jbd,il}, flp] = ...
                    hssupd(lh,V{wh},rk{wh},SI{jbd,jbd,c2},-SI{jbd,ec2,c2}*tmp{2,2},sh,ih);
                    
                    flp2 = flp2 + flp + 2*ncol2*(numel(C{2,2})+numel(tmp{2,2}));
                end
                
                for ibd = bds
                    if(ibd==jbd); continue; end
                    nrow1 = min( rk{whb(ibd,c1tr)}(shb(ibd,c1tr),ihb(ibd,c1tr)), szc(1,ibd) );
                    nrow2 = min( rk{whb(ibd,c2tr)}(shb(ibd,c2tr),ihb(ibd,c2tr)), szc(2,ibd) );
                    SI{ibd,jbd,il} = zeros(nrow1+nrow2,ncol1+ncol2);
                    % c1, c1
                    if (nrow1*ncol1 > 0)
                        SI{ibd,jbd,il}(1:nrow1,1:ncol1) = ...
                        SI{ibd,jbd,c1} - SI{ibd,ec1,c1} * tmp{1,1};
                        flp2 = flp2 + 2*nrow1*ncol1*size(tmp{1,1},1);
                    end
                    % c1, c2
                    if (nrow1*ncol2 > 0)
                        SI{ibd,jbd,il}(1:nrow1,ncol1+1:ncol1+ncol2) = ...
                        -SI{ibd,ec1,c1} * tmp{1,2};

                        flp2 = flp2 + 2*nrow1*ncol2*size(tmp{1,2},1);
                    end
                    % c2, c1
                    if (nrow2*ncol1 > 0)
                        SI{ibd,jbd,il}(nrow1+1:nrow1+nrow2,1:ncol1) = ...
                        -SI{ibd,ec2,c2} * tmp{2,1};
                        flp2 = flp2 + 2*nrow2*ncol1*size(tmp{2,1},1);
                    end
                    
                    % c2, c2
                    if (nrow2*ncol2 > 0)
                        SI{ibd,jbd,il}(nrow1+1:nrow1+nrow2,ncol1+1:ncol1+ncol2) = ...
                        SI{ibd,jbd,c2} - SI{ibd,ec2,c2} * tmp{2,2};
                        flp2 = flp2 + 2*nrow2*ncol2*size(tmp{2,2},1);
                    end
                end
                
                % 1 level of HSS column compression, skip this at root
                totsz = sum(sz(:,il));
                if ( lhb(jbd,iltr) < lhb(jbd,c1tr) && totsz > sz(jbd,il))
                    [U1,sig,V1]=svd(vertcat(SI{[1:jbd-1,jbd+1:nbd],jbd,il}));
                    sig = diag(sig);
                    flp2 = flp2 + 2*numel(V1)*size(U1,1);
                    if (sig(end)<=sig(1)*tol)
                        % rank deficient
                        rkt = max(find(sig<=sig(1)*tol,1)-1, 1);
                    else
                        rkt = length(sig);
                    end
                    rk{wh}(sh,ih) = rkt;
                    V{wh}{2,sh,ih} = V1;
                    stg2 = stg2 + size(V1,1)*rkt;
                    for ibd = bds
                        if(ibd==jbd); continue; end
                        SI{ibd,jbd,il} = SI{ibd,jbd,il}*V1(:,1:rkt);
                    end
                end
                
            end
            
            for ibd = bds
                if ( lhb(ibd,iltr) < lhb(ibd,c1tr) && totsz > sz(ibd,il))
                    [U1,sig,V1]=svd(horzcat(SI{ibd,[1:ibd-1,ibd+1:nbd],il}));
                    sig = diag(sig);
                    wh = whb(ibd,iltr);
                    sh = shb(ibd,iltr);
                    ih = ihb(ibd,iltr);
                    rkt = rk{wh}(sh,ih);
                    V{wh}{1,sh,ih} = U1;
                    stg2 = stg2 + size(U1,1)*rkt;
                    for jbd = bds
                        if(ibd==jbd); continue; end
                        SI{ibd,jbd,il} = U1(:,1:rkt)'*SI{ibd,jbd,il};
                    end
                end
                stg2 = stg2 + numel(SI{ibd,ec1,c1}) + numel(SI{ec1,ibd,c1})...
                            + numel(SI{ibd,ec2,c2}) + numel(SI{ec2,ibd,c2});
            end
        end
    end
end

t2 = toc;

counts = struct('time',[t1,t2],'nflops',[flp1,flp2],'storage',[stg1,stg2]);
end

function [Bo, flp] = hssupdint2(nlvl,V,rk,Bi1,Bi2,M,sh,ih)
n1d   = 2^nlvl;
Bo    = cell(1, n1d*3-1);
Bo{2} = M{1,2};
Bo{3} = M{2,1};
if (nlvl == 1)
    % diagonal
    Bo{4} = Bi1 + V{1,sh,ih*2  }(:,1:size(M{1,1},1)) * M{1,1} * V{2,sh,ih*2  }(:,1:size(M{1,1},2))';
    Bo{5} = Bi2 + V{1,sh,ih*2+1}(:,1:size(M{2,2},1)) * M{2,2} * V{2,sh,ih*2+1}(:,1:size(M{2,2},2))';
    flp = 2*(size(Bi1,1)*numel(M{1,1}) + numel(Bi1)*size(M{1,1},2)...
            +size(Bi2,1)*numel(M{2,2}) + numel(Bi2)*size(M{2,2},2));
    return;
end
Bo{n1d*2}   = M{1,1};
Bo{n1d*2+1} = M{2,2};
flp = 0;

nod  = 2;
n1v  = ih*2;
for lvl = 2:nlvl
    for ib = nod*2-1:-1:nod
        %Bo lvl+1 <- Bi lvl
        Bo{ib+nod}   = Bi1{ib};
        Bo{ib+nod*2} = Bi2{ib};
    end
    
    for ib = nod*2-1:-1:nod
        id  = ib+n1d*2-nod;
        iv  = ib+n1v-nod;
        
        T  = V{1,sh,iv}(:,1:size(Bo{id},1))*Bo{id}*(V{2,sh,iv}(:,1:size(Bo{id},2)))';
        flp = flp + 2*(size(T,1)*numel(Bo{id}) + numel(T)*size(Bo{id},2));
        
        c1b = ib*2;
        c2b = ib*2+1;
        c1d = (ib-nod+n1d)*2;
        c2d = c1d+1;
        
        rk1 = rk(sh,iv*2);
        rk2 = rk1;
        Bo{c2b} = Bo{c2b} + T(rk2+1:end,1:rk1);
        Bo{c1b} = Bo{c1b} + T(1:rk1,rk2+1:end);
        Bo{c2d} = T(rk1+1:end,rk2+1:end);
        Bo{c1d} = T(1:rk1,1:rk2);
    end
    nod = nod*2;
    n1v = n1v*2;
end

for ib = 0:nod-1
    id  = ib+n1d*2;
    iv  = ib+n1v;
    
    T = V{1,sh,iv}(:,1:size(Bo{id},1))*Bo{id}*(V{2,sh,iv}(:,1:size(Bo{id},2)))';
    flp = flp + 2*(size(T,1)*numel(Bo{id}) + numel(T)*size(Bo{id},2));
    
    if (id-n1d<=length(Bi1))
        Bo{id} = T + Bi1{id-n1d};
    else
        Bo{id} = T + Bi2{id-length(Bi1)-1};
    end
end

end