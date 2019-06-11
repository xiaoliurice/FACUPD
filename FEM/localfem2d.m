function [kit, A, M, Mb, pd, A0] = localfem2d(P,Pb,Bvtx,hvtx,bvps,bghs, nlvl,nsub,lvb,clf, kloc, hloc)
% local problem A*u = M*f+Mb*g

% Build reference element
[x,y] = Nodes2D(P);
[r,s] = xytors(x,y);

V2 = Vandermonde2D(P,r,s);
[Dr,Ds] = Dmatrices2D(P,r,s,V2);
mass2 = inv(V2);
mass2 = mass2'*mass2;

V1 = Vandermonde1D(P, r(1:P+1));
inV1 = inv(V1);
mass1 = inV1'*inV1;

% Data structures
% --from vertex (iv) to coordinates (x,y) and nodal id
% Xvtx(iv), Yvtx(iv), v2n(iv)
% --from face (if) to vertices, nodals, Jacobians
% f2v(if,1:2), n4f(:,if), detJ4f(if)
% --from element (ie) to vertices, nodals, Jacobians
% e2v(ie,1:3), n4e(:,ie), invJ4e(:,:,ie), detJ4e(ie)
% --from nodal (in) to coordinate (x,y)
% Xndl(in), Yndl(in)

% size etc
Bndl = Bvtx*P;
Nvtx = prod(Bvtx+1); Nndl = prod(Bndl+1);
P1 = P+1; P2 = P1*(P1+1)/2; Pb1 = Pb+1;
nzfac  = P1*P1; nzfacb = P1*Pb1; nzele  = P2*P2;
bd   = [Bvtx(2), Bvtx(2), Bvtx(1), Bvtx(1)];
nbd  = length(bd);
Nfac = sum(bd);

% Xvtx(Nx,Ny), Yvtx(Nx,Ny), vtx2ndl(Nx,Ny) all integers, ele2vtx(Nele,3), ndl4ele(P2,Nele)
[Yvtx, Xvtx] = meshgrid(0:Bvtx(2), 0:Bvtx(1));
v2n = (Yvtx*(Bndl(1)+1)+Xvtx)*P+1;
e2v = delaunay(Xvtx*hvtx(1),Yvtx*hvtx(2));
Nele = size(e2v,1);

% stiffness and mass matrix
row = zeros(nzele*Nele+nzfac*Nfac,1);
col = zeros(nzele*Nele+nzfac*Nfac,1);
stf = zeros(nzele*Nele,1);
mas = zeros(nzele*Nele+nzfac*Nfac,1);

Xndl = zeros(Bndl'+1);
Yndl = zeros(Bndl'+1);
ne1 = zeros(P2,1);
ne2 = zeros(P2,1);
ee = 1;
for j = 0:P
    for i = 0:P-j
        ne1(ee) = i;
        ne2(ee) = j;
        ee = ee + 1;
    end
end

oe = 0;
for ee = 1:Nele
    oee = oe + nzele;
    
    c0 = v2n(e2v(ee,1));
    c1 = (v2n(e2v(ee,2))-c0)/P;
    c2 = (v2n(e2v(ee,3))-c0)/P;
    n4e = c1*ne1 + c2*ne2 + c0;
    
    x = 0.5*(-(r+s)*Xvtx(e2v(ee,1))+(1+r)*Xvtx(e2v(ee,2))+(1+s)*Xvtx(e2v(ee,3)))*hvtx(1);
    y = 0.5*(-(r+s)*Yvtx(e2v(ee,1))+(1+r)*Yvtx(e2v(ee,2))+(1+s)*Yvtx(e2v(ee,3)))*hvtx(2);
    Xndl(n4e) = x;
    Yndl(n4e) = y;
    row(oe+1:oee) = repmat(n4e,P2,1);
    col(oe+1:oee) = kron(n4e,ones(P2,1));
    
    invJ = [Dr(1,:)*x,Ds(1,:)*x;Dr(1,:)*y,Ds(1,:)*y]; 
    detJ = det(invJ);  % |d(x,y)/d(r,s)|
    invJ = inv(invJ);  %  d(r,s)/d(x,y)
    Dx = Dr*invJ(1,1)+Ds*invJ(2,1);
    Dy = Dr*invJ(1,2)+Ds*invJ(2,2);
    stf(oe+1:oee) = reshape(detJ*(Dx'*mass2*Dx + Dy'*mass2*Dy), [nzele,1]);
    mas(oe+1:oee) = detJ*mass2(:);
    
    oe = oee;
end

% fac2vtx(Nfac,2)
f2v = horzcat([(0:Bvtx(2)-1)*(Bvtx(1)+1)+1,...
               (1:Bvtx(2)  )*(Bvtx(1)+1), ...
                1:Bvtx(1), ...
               Nvtx-Bvtx(1):Nvtx-1]',...
               [(1:Bvtx(2)  )*(Bvtx(1)+1)+1,...
                (2:Bvtx(2)+1)*(Bvtx(1)+1), ...
                 2:Bvtx(1)+1, ...
               Nvtx-Bvtx(1)+1:Nvtx]');

rowb = zeros(nzfacb*Nfac,1);
colb = zeros(nzfacb*Nfac,1);
masb = zeros(nzfacb*Nfac,1);

of = 0;
for ee = 1:Nfac
    oee = oe + nzfac;
    off = of + nzfacb;
    
    c0 = v2n(f2v(ee,1));
    c1 = round( (v2n(f2v(ee,2))-c0) / P );
    n4f = (c0:c1:c0+P*c1)';
    
    i = v2n(f2v(ee,1));
    j = v2n(f2v(ee,2));
    detJ = norm([Xndl(i)-Xndl(j),Yndl(i)-Yndl(j)])/2; % 1D Jacobian d|xj-xi|/d|1-(-1)|\
    
    row(oe+1:oee) = repmat(n4f,P1,1);
    col(oe+1:oee) = kron(n4f,ones(P1,1));
    mas(oe+1:oee) = detJ*mass1(:);
    
    rowb(of+1:off) = repmat(n4f,Pb1,1);
    colb(of+1:off) = kron((ee-1)*Pb1+1:ee*Pb1,ones(P1,1));
    masb(of+1:off) = sqrt(detJ)*reshape(inV1(1:Pb1,:)',nzfacb,1);
    
    oe = oee;
    of = off;
end

M = sparse(row(1:nzele*Nele), col(1:nzele*Nele), mas(1:nzele*Nele));
Mb = cell(nbd,1);
ob = 0; of = 0;
for ibd = 1:nbd
    obb = ob + bd(ibd)*Pb1; off = of + bd(ibd)*nzfacb;
    Mb{ibd} = sparse(rowb(of+1:off),colb(of+1:off)-ob,masb(of+1:off),Nndl,obb-ob);
    ob = obb; of = off;
end

nod  = 2^(nlvl-1);
A    = cell(1,nod);
A0   = cell(1,nod);
pd   = zeros(Nndl,nod);

na = Bvtx.*P;
kit  = struct('x',zeros((na.*nsub+1)'), 'y',zeros((na.*nsub+1)'), 'k',zeros((na.*nsub+1)'));
kiti = reshape( 1:length(kit.x(:)), (na.*nsub+1)' );

ktmp = zeros(size(Xndl));
for ii = 1:nod
    il = ii-1+nod;
    bvpi = [2, 2, 2, 2]; % Robin

    for ibd = 1:nbd
        if (lvb(ibd,il) == 1)
            % global boundary
            bvpi(ibd) = bvps(ibd);
        end
    end

    % sparse matrices
    aval = zeros(size(mas));
    oe = 0;
    for ee = 1:Nele
        c0 = v2n(e2v(ee,1));
        c1 = (v2n(e2v(ee,2))-c0)/P;
        c2 = (v2n(e2v(ee,3))-c0)/P;
        n4e = c1*ne1 + c2*ne2 + c0;
        
        % interp nodals from a slice of the model
        x = max(1, min( 1+(Xndl(n4e)+kloc{ii}.d(1))/hloc(1), size(kloc{ii}.val,1) ));
        y = max(1, min( 1+(Yndl(n4e)+kloc{ii}.d(2))/hloc(2), size(kloc{ii}.val,2) ));
        kln = interp2(kloc{ii}.val,y,x);
        oee = oe + nzele;
        ktmp(n4e) = kln;
        aval(oe+1:oee) = stf(oe+1:oee) - mas(oe+1:oee) .* kron(kln, kln);
        oe = oee;
    end

    for ibd = 1:nbd
        bb  = bvpi(ibd);
        tmp = bghs(1,1,bb)/bghs(1,2,bb);
        oee = oe + bd(ibd)*nzfac;
        aval(oe+1:oee) = tmp*mas(oe+1:oee);
        oe = oee;
    end
    A{ii} = sparse(row,col,aval);
    
    % remove artificial boundary terms to get the original matrix
    oe = nzele*Nele;
    for ibd = 1:nbd
        oee = oe + bd(ibd)*nzfac;
        if lvb(ibd,il) ~= 1
            aval(oe+1:oee) = 0;
        end
        oe = oee;
    end
    A0{ii} = sparse(row,col,aval);
    
    fst1 = (clf(1,ii)-1)*na(1)+1;
    lst1 =  clf(1,ii)   *na(1)+1;
    fst2 = (clf(2,ii)-1)*na(2)+1;
    lst2 =  clf(2,ii)   *na(2)+1;
    p = kiti(fst1:lst1,fst2:lst2);
    pd(:,ii) = p(:);
    kit.x(p) = Xndl + (clf(1,ii)-1)* (Bvtx(1)*hvtx(1));
    kit.y(p) = Yndl + (clf(2,ii)-1)* (Bvtx(2)*hvtx(2));
    kit.k(p) = ktmp;
end
end