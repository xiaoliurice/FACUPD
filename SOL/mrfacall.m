function [FI,SI,FE,SE,sz, counts] = mrfacall(nlvl,lvb,bvps,bghs, A,Mb)
% construct every interior and exterior factor

% interior
[FI,SI,sz, counts,c1,c2,ec1,ec2] = mrfacint(1,nlvl,lvb,bvps,bghs, A,Mb);

nbd = 4;
nod = 2^nlvl-1;

% exterior
FE = cell(nbd,nbd+2,nod);
SE = cell(nbd,nbd,nod);

flp2 = 0.0;
stg2 = 0.0;

SE{ec1,ec1,c1} = SI{ec2,ec2,c2};
SE{ec2,ec2,c2} = SI{ec1,ec1,c1};

tic;
for lvl = 2:nlvl-1
    nod = 2^(lvl-1);
    for il = nod:nod*2-1
        c1 = il*2;
        c2 = il*2+1;
        
        % size for the assembly
        szc = zeros(2,nbd);
        for ibd = 1:nbd
            szc(:,ibd) = [sz(ibd,c1),sz(ibd,c2)];
            if (lvb(ibd,c1) == lvl+1)
                szc(1,ibd) = 0;
                ec1 = ibd;
            end
            if (lvb(ibd,c2) == lvl+1)
                szc(2,ibd) = 0;
                ec2 = ibd;
            end
        end
        bds1 = find(szc(1,:));
        bds2 = find(szc(2,:));
        bds1 = bds1(:)';
        bds2 = bds2(:)';

        % neighbor of c2
        % N: SE{:,:,il}, M: SI{:,:,c2}
        % NM-I;
        ne = sum(szc(2,bds2));
        for jbd = bds2
            for ibd = bds2
                if (ibd==jbd)
                   FE{jbd,ibd,c2} = -eye(szc(2,jbd));
                else
                   FE{jbd,ibd,c2} = zeros(szc(2,jbd),szc(2,ibd));
                end
                for kbd = bds2
                   FE{jbd,ibd,c2} = FE{jbd,ibd,c2} ...
                  +SE{jbd,kbd,il}(szc(1,jbd)+1:end,szc(1,kbd)+1:end)*SI{kbd,ibd,c2};
                end
            end
        end


        % block LU
        for j = 1:length(bds2)
            jbd = bds2(j);
            % LU of A11
            [FE{jbd,nbd+1,c2},FE{jbd,nbd+2,c2}] = lu(FE{jbd,jbd,c2});
            FE{jbd,jbd,c2} = [];
            for i = j+1:length(bds2)
                ibd = bds2(i);
                % inv(L)*A12
                FE{jbd,ibd,c2} = FE{jbd,nbd+1,c2} \ FE{jbd,ibd,c2};
                % A21*inv(U)
                FE{ibd,jbd,c2} = FE{ibd,jbd,c2} / FE{jbd,nbd+2,c2};
            end
            for i = j+1:length(bds2)
                ibd = bds2(i);
                for k = j+1:length(bds2)
                   kbd = bds2(k);
                   FE{kbd,ibd,c2} = FE{kbd,ibd,c2} - FE{kbd,jbd,c2}*FE{jbd,ibd,c2};
                end
            end
        end
        flp2 = flp2 + (2+2/3)*ne^3;
        stg2 = stg2 + ne^2;
        
        % / KN  \
        % \MKN-I/
        tmp = cell(nbd,2);
        for ibd = bds2
            tmp{ibd,1} = zeros(szc(2,ibd),sz(ec2,c2));
            for kbd = bds2
                tmp{ibd,1} = tmp{ibd,1} ...
               +SE{ibd,kbd,il}(szc(1,ibd)+1:end,szc(1,kbd)+1:end)*SI{kbd,ec2,c2};
            end
        end

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

        for ibd = bds2
            tmp{ibd,2} = -SI{ibd,ec2,c2};
            for kbd = bds2
                tmp{ibd,2} = tmp{ibd,2} + SI{ibd,kbd,c2}*tmp{kbd,1};
            end
        end
        flp2 = flp2 + sz(ec2,c2)*(6*ne^2);

        % interface c1
        SE{ec1,ec1,c1} = SI{ec2,ec2,c2};
        for kbd = bds2
            SE{ec1,ec1,c1} = SE{ec1,ec1,c1} - SI{ec2,kbd,c2}*tmp{kbd,1};
        end
        stg2 = stg2 + sz(ec1,c1)^2;
        flp2 = flp2 + 2*ne*sz(ec1,c1)^2;
        
        % neigbor of c1, interface of c1
        for ibd = bds1
            SE{ibd,ec1,c1} = zeros(szc(1,ibd),sz(ec1,c1));
            for kbd = bds2
                SE{ibd,ec1,c1} = SE{ibd,ec1,c1}...
              - SE{ibd,kbd,il}(1:szc(1,ibd),szc(1,kbd)+1:end)*tmp{kbd,2};
            end
            stg2 = stg2 + szc(1,ibd)*sz(ec1,c1);
            flp2 = flp2 + 2*ne*szc(1,ibd)*sz(ec1,c1);
        end

        for jbd = bds1
            tmp = cell(nbd,2);
            % / K\
            % \MK/
            for ibd = bds2
                tmp{ibd,1} = SE{ibd,jbd,il}(szc(1,ibd)+1:end,1:szc(1,jbd));
            end 
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
            for ibd = bds2
                tmp{ibd,2} = zeros(szc(2,ibd),szc(1,jbd));
                for kbd = bds2
                    tmp{ibd,2} = tmp{ibd,2}+SI{ibd,kbd,c2}*tmp{kbd,1};
                end
            end
            flp2 = flp2 + szc(1,jbd)*(4*ne^2);

            % interface c1, neighbor c1
            SE{ec1,jbd,c1} = zeros(sz(ec1,c1),szc(1,jbd));
            for kbd = bds2
                SE{ec1,jbd,c1} = SE{ec1,jbd,c1} - SI{ec2,kbd,c2}*tmp{kbd,1};
            end
            stg2 = stg2 + sz(ec1,c1)*szc(1,jbd);
            flp2 = flp2 + 2*ne*sz(ec1,c1)*szc(1,jbd);
            

            % neigbor c1, neigbor c1
            for ibd = bds1
                SE{ibd,jbd,c1} = SE{ibd,jbd,il}(1:szc(1,ibd),1:szc(1,jbd));
                for kbd = bds2
                    SE{ibd,jbd,c1} = SE{ibd,jbd,c1}...
                   -SE{ibd,kbd,il}(1:szc(1,ibd),szc(1,kbd)+1:end) * tmp{kbd,2};
                end
                stg2 = stg2 + szc(1,ibd)*szc(1,jbd);
                flp2 = flp2 + 2*ne*szc(1,ibd)*szc(1,jbd);
            end
           
        end

        % neighbor of c1
        % N: SE{il}, M: SI{c1}
        % NM-I;
        ne = sum(szc(1,bds1));
        for jbd = bds1
            for ibd = bds1
                if (ibd==jbd) 
                   FE{jbd,ibd,c1} = -eye(szc(1,jbd));
                else
                   FE{jbd,ibd,c1} = zeros(szc(1,jbd),szc(1,ibd));
                end
                for kbd = bds1
                   FE{jbd,ibd,c1} = FE{jbd,ibd,c1} ...
                  +SE{jbd,kbd,il}(1:szc(1,jbd),1:szc(1,kbd))*SI{kbd,ibd,c1};
                end
            end
        end

        % block LU
        for j = 1:length(bds1)
            jbd = bds1(j);
            % LU of A11
            [FE{jbd,nbd+1,c1},FE{jbd,nbd+2,c1}] = lu(FE{jbd,jbd,c1});
            FE{jbd,jbd,c1} = [];
            for i = j+1:length(bds1)
                ibd = bds1(i);
                % inv(L)*A12
                FE{jbd,ibd,c1} = FE{jbd,nbd+1,c1} \ FE{jbd,ibd,c1};
                % A21*inv(U)
                FE{ibd,jbd,c1} = FE{ibd,jbd,c1}   / FE{jbd,nbd+2,c1};
            end
            for i = j+1:length(bds1)
                ibd = bds1(i);
                for k = j+1:length(bds1)
                   kbd = bds1(k);
                   FE{kbd,ibd,c1} = FE{kbd,ibd,c1} - FE{kbd,jbd,c1}*FE{jbd,ibd,c1};
                end
            end
        end
        flp2 = flp2 + (2+2/3)*ne^3;
        stg2 = stg2 + ne^2;
        
        % / KN  \
        % \MKN-I/
        tmp = cell(nbd,2);
        for ibd = bds1
            tmp{ibd,1} = zeros(szc(1,ibd),sz(ec1,c1));
            for kbd = bds1
                tmp{ibd,1} = tmp{ibd,1} ...
               +SE{ibd,kbd,il}(1:szc(1,ibd),1:szc(1,kbd))*SI{kbd,ec1,c1};
            end
        end

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

        for ibd = bds1
            tmp{ibd,2} = -SI{ibd,ec1,c1};
            for kbd = bds1
                tmp{ibd,2} = tmp{ibd,2} + SI{ibd,kbd,c1}*tmp{kbd,1};
            end
        end
        flp2 = flp2 + sz(ec1,c1)*(6*ne^2);

        % interface c2
        SE{ec2,ec2,c2} = SI{ec1,ec1,c1};
        for kbd = bds1
            SE{ec2,ec2,c2} = SE{ec2,ec2,c2} - SI{ec1,kbd,c1}*tmp{kbd,1};
        end
        stg2 = stg2 + sz(ec2,c2)^2;
        flp2 = flp2 + 2*ne*sz(ec2,c2)^2;
        
        % neigbor of c2, interface of c2
        for ibd = bds2
            SE{ibd,ec2,c2} = zeros(szc(2,ibd),sz(ec2,c2));
            for kbd = bds1
                SE{ibd,ec2,c2} = SE{ibd,ec2,c2}...
               -SE{ibd,kbd,il}(szc(1,ibd)+1:end,1:szc(1,kbd))*tmp{kbd,2};
            end
            stg2 = stg2 + szc(2,ibd)*sz(ec2,c2);
            flp2 = flp2 + 2*ne*szc(2,ibd)*sz(ec2,c2);
        end

        for jbd = bds2
            tmp = cell(nbd,2);
            % / K\
            % \MK/
            for ibd = bds1
                tmp{ibd,1} = SE{ibd,jbd,il}(1:szc(1,ibd),szc(1,jbd)+1:end);
            end 
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
            for ibd = bds1
                tmp{ibd,2} = zeros(szc(1,ibd),szc(2,jbd));
                for kbd = bds1
                    tmp{ibd,2} = tmp{ibd,2}+SI{ibd,kbd,c1}*tmp{kbd,1};
                end
            end
            flp2 = flp2 + szc(2,jbd)*(4*ne^2);

            % interface c2, neighbor c2
            SE{ec2,jbd,c2} = zeros(sz(ec2,c2),szc(2,jbd));
            for kbd = bds1
                SE{ec2,jbd,c2} = SE{ec2,jbd,c2} - SI{ec1,kbd,c1}*tmp{kbd,1};
            end
            stg2 = stg2 + sz(ec2,c2)*szc(2,jbd);
            flp2 = flp2 + 2*ne*sz(ec2,c2)*szc(2,jbd);
            

            % neigbor c2, neigbor c2
            for ibd = bds2
                SE{ibd,jbd,c2} = SE{ibd,jbd,il}(szc(1,ibd)+1:end,szc(1,jbd)+1:end);
                for kbd = bds1
                    SE{ibd,jbd,c2} = SE{ibd,jbd,c2}...
                   -SE{ibd,kbd,il}(szc(1,ibd)+1:end,1:szc(1,kbd)) * tmp{kbd,2};
                end
                stg2 = stg2 + szc(2,ibd)*szc(2,jbd);
                flp2 = flp2 + 2*ne*szc(2,ibd)*szc(2,jbd);
            end
           
        end

    end
end
t2=toc;

counts = struct('time',[sum(counts.time), t2],'nflops',[sum(counts.nflops),flp2],'storage',[sum(counts.storage),stg2]);
end