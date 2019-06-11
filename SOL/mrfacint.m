function [FI,SI,sz, counts, c1,c2, ec1,ec2] = mrfacint(rid,nlvl,lvb,bvps,bghs, A,Mb)
% construct interior factors
% SI: Schur complement as interior boundary maps
% FI: LU factors
% rid: ID of root  (1,1 for root)

rlvl = floor(log2(rid)) + 1;

nbd = 4;
nod = 2^(nlvl-rlvl+1)-1; % number of nodes in the subtree
sz = zeros(nbd,nod);  % size of each BD components

% interior
FI = cell(4,nod); % L,U,p,q
SI = cell(nbd,nbd,nod);


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
    stg2 = stg2 + sum(sz(:,il))^2;
end
t1=toc;

%%%%%%%%%%%%% merge uncompressed levels %%%%%%%%%%%%%%%%%
tic;
for lvl = nlvl-1:-1:rlvl
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
        stg2 = stg2 + ne^2 + sum(sz(:,il))^2; %F, S
        flp2 = flp2 + (2/3+2)*ne^3;
        
        bds = find(sz(:,il)); % nonzero blocks
        bds = bds(:)';
        for jbd = bds
            tmp = cell(2,2);
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


t2 = toc;

counts = struct('time',[t1,t2],'nflops',[flp1,flp2],'storage',[stg1,stg2]);
end