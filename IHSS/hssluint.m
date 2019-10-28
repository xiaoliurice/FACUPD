function [L,U, L2,U2, C, flp,stg ] = hssluint(nlvl,V,rk,S1,S2)
flp = 0.0;
stg = 0.0;
C = cell(2,2);
if(nlvl==1)
    ih = 1;
    % only outer level
    ne = size(S1,1);
    stg = stg + ne^2;
    flp = flp + (2/3+2)*ne^3;
    [L, U] = lu(S2*S1 - eye(ne));
    if (rk(1,ih)>0)
        % inv2x2 * [Uc1;0]
        C{1,1} = U \ (L \ (S2*V{1,1,ih}(:,1:rk(1,ih))));
        if(rk(2,ih)>0)
            C{2,1} = (V{2,2,ih}(:,1:rk(2,ih)))'*(S1*C{1,1}-V{1,1,ih}(:,1:rk(1,ih)));
            flp = flp + rk(1,ih)*(2*ne^2 + 2*ne*rk(2,ih));
        end
        C{1,1} = (V{2,1,ih}(:,1:rk(1,ih)))'*C{1,1};
        flp = flp + rk(1,ih)*(4*ne^2 + 2*ne*rk(1,ih));
    end
    if (rk(2,ih)>0)
        % inv2x2 * [0;Uc2]
        C{1,2} = U \ (L \ V{1,2,ih}(:,1:rk(2,ih)));
        C{2,2} = (V{2,2,ih}(:,1:rk(2,ih)))'*S1*C{1,2};
        flp = flp + rk(2,ih)*(4*ne^2+2*ne*rk(2,ih));
        if (rk(1,ih)>0)
            C{1,2} = (V{2,1,ih}(:,1:rk(1,ih)))'*C{1,2};
            flp = flp + rk(2,ih)*(2*ne*rk(1,ih));
        end
    end
    L2 = S1;
    U2 = S2;
else
    % inner & outer levels
    nod = 2^(nlvl-1);
    n1d = 2^nlvl;
    L  = cell(1,n1d-1);
    U  = cell(1,n1d-1);
    L2 = cell(1,n1d-1);
    U2 = cell(1,n1d-1);
    D  = cell(2,2,nod);
    for lvl=nlvl:-1:1
        for iv = nod:nod*2-1
            it = iv+1-nod;
            if(lvl==nlvl)
                id = iv+n1d-nod;
                
                D{1,1,it} =  V{1,1,iv}'*S1{id}*V{2,1,iv};
                D{2,2,it} =  V{1,2,iv}'*S2{id}*V{2,2,iv};
                D{1,2,it} = -V{1,1,iv}'*V{2,2,iv};
                D{2,1,it} = -V{1,2,iv}'*V{2,1,iv};
                szt = [size(S1{id},1),size(S2{id},1)];
                
                flp = flp + 2 * szt(1)^3 * 6;
            else
                c1t = it*2-1;
                c2t = it*2;
                c1v = iv*2;
                c2v = iv*2+1;
                
                D{1,1,it} = [D{1,1,c1t},S1{c1v};S1{c2v},D{1,1,c2t}];
                D{2,2,it} = [D{2,2,c1t},S2{c1v};S2{c2v},D{2,2,c2t}];
                D{1,2,it} = blkdiag(D{1,2,c1t},D{1,2,c2t});
                D{2,1,it} = blkdiag(D{2,1,c1t},D{2,1,c2t});
                szt = [size(D{1,1,it},1),size(D{2,2,it},1)];       
                
                if(lvl>1)
                    D{1,1,it} = V{1,1,iv}'*D{1,1,it}*V{2,1,iv};
                    D{2,2,it} = V{1,2,iv}'*D{2,2,it}*V{2,2,iv};
                    D{1,2,it} = V{1,1,iv}'*D{1,2,it}*V{2,2,iv};
                    D{2,1,it} = V{1,2,iv}'*D{2,1,it}*V{2,1,iv};
                
                    flp = flp + 4*sum(szt.^2)*sum(szt);
                end
            end
            
            if(lvl>1)
                
                rkt = rk(:,iv);
                % leading block A11 = L1*U1
                [L{iv},U{iv}] = lu(...
                  [D{1,1,it}(rkt(1)+1:end,rkt(1)+1:end),D{1,2,it}(rkt(1)+1:end,rkt(2)+1:end);...
                   D{2,1,it}(rkt(2)+1:end,rkt(1)+1:end),D{2,2,it}(rkt(2)+1:end,rkt(2)+1:end)]);
                % U1\A21
                L2{iv} = ...
                  [D{1,1,it}(1:rkt(1),rkt(1)+1:end),D{1,2,it}(1:rkt(1),rkt(2)+1:end);...
                   D{2,1,it}(1:rkt(2),rkt(1)+1:end),D{2,2,it}(1:rkt(2),rkt(2)+1:end)]/U{iv};
                % L1\A12
                U2{iv} = L{iv}\...
                  [D{1,1,it}(rkt(1)+1:end,1:rkt(1)),D{1,2,it}(rkt(1)+1:end,1:rkt(2));...
                   D{2,1,it}(rkt(2)+1:end,1:rkt(1)),D{2,2,it}(rkt(2)+1:end,1:rkt(2))];
                n2 = sum(rkt);
                n1 = sum(szt)-n2;
                
                % Schur complement
                D{1,1,it} = D{1,1,it}(1:rkt(1),1:rkt(1)) ...
                  - L2{iv}(1:rkt(1),:)*U2{iv}(:,1:rkt(1));
                D{1,2,it} = D{1,2,it}(1:rkt(1),1:rkt(2)) ...
                  - L2{iv}(1:rkt(1),:)*U2{iv}(:,rkt(1)+1:end);
                D{2,1,it} = D{2,1,it}(1:rkt(2),1:rkt(1)) ...
                  - L2{iv}(rkt(1)+1:end,:)*U2{iv}(:,1:rkt(1));
                D{2,2,it} = D{2,2,it}(1:rkt(2),1:rkt(2)) ...
                  - L2{iv}(rkt(1)+1:end,:)*U2{iv}(:,rkt(1)+1:end);
                flp = flp + 2/3*n1^3 + 2*n1^2*n2 + 2*n2^2*n1;
                stg = stg + n1^2 + 2*n1*n2;
            else
                n1 = sum(szt);
                [L{iv},U{iv}] = lu([D{1,1,it},D{1,2,it};D{2,1,it},D{2,2,it}]);
                flp = flp + 2/3*n1^3;
                stg = stg + n1^2;
                if(rk(1,iv)>0)
                    T1 = U{iv}\ (L{iv} \...
                        [V{1,1,iv}(:,1:rk(1,iv));zeros(szt(2),rk(1,iv))] );
                    C{1,1} = V{2,1,iv}(:,1:rk(1,iv))'*T1(1:szt(1),:);
                end
                
                if(rk(2,iv)>0)
                    %2
                    T2 = U{iv}\ (L{iv} \...
                        [zeros(szt(1),rk(2,iv));V{1,2,iv}(:,1:rk(2,iv) )]);
                    C{2,2} = V{2,2,iv}(:,1:rk(2,iv))'*T2(szt(1)+1:end,:);
                end
                flp = flp + 2*n1^2*sum(rk(:,iv));
                
                if(rk(1,iv)*rk(2,iv)>0)
                    %1,2
                    C{2,1} = V{2,2,iv}(:,1:rk(2,iv))'*T1(szt(1)+1:end,:);
                    C{1,2} = V{2,1,iv}(:,1:rk(1,iv))'*T2(1:szt(1),:);
                    flp = flp + numel(C{2,1})*szt(2) + numel(C{1,2})*szt(1);
                end
                
            end
        end
        nod = nod/2;
    end
end
end