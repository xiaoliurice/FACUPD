%% main script
addpath('FEM','SOL');

%% setup the domain
B0=[80;80];   % num of vertices
nlvl=7;       % num of levels of domain partitioning
[nsub,lvb,clf] = bisect(B0,nlvl);    
[B0, Bvtx] = partition(nsub, 0, B0);
Pele = 4; Pfac = 2;  % polynomial order of interior and boundary basis functions
Bndl = B0*Pele+1; Bmul = Bvtx*Pele; hvtx = [Pfac;Pfac];
fprintf('--Setup--\n  Number of vertices: global %dx%d, local %dx%d\n',B0,Bvtx);
fprintf('  Number of nodal points: %dx%d\n', Bndl);

% setup the test problem
kmax = 2*pi/8;             % dimensionless wavenumber: 2pi / (sampling rate)
kloc = cell(1,prod(nsub));
for i = 1:size(clf,2)
    kloc{i} = struct('val',ones(Bvtx'+1)*kmax,'d',zeros(size(Bvtx)));
end

% boundary condition: Up Down Left Right
bvps =[2,2,2,2];          % 1 Neumann-to-Dirichlet; 2 Impedance-to-Impedance
bghs = zeros(2,2,2);      % coefficients
bghs(:,:,1) = [0,1;1,0];
bghs(:,:,2) = [-1i*kmax,1; -1i*kmax,-1];

% generate the matrix
tic;
[kit,A1, M,Mb,pd] = localfem2d(Pele,Pfac,Bvtx,hvtx,bvps,bghs, nlvl,nsub,lvb,clf, kloc,hvtx);
t=toc;
fprintf('  Matrix: n %d, nnz %d; FEM time %fs\n',prod(Bndl), nnz(A1{1})*length(A1),t);

disp('--Factorization of background--');
[FI,SI,FE,SE,sz, counts] = mrfacall(nlvl,lvb,bvps,bghs, A1,Mb);
fprintf('  Total time: %fs, nflops: %.2e, storage %.2e\n', sum(counts.time), sum(counts.nflops), sum(counts.storage));

% right hand side
Bsrc = B0.*hvtx.*[0.5;0.5];
r = zeros(size(M,1), prod(nsub));
for i = 1:size(r,2)
    r(:,i) = M*exp( - (kit.x(pd(:,i))-Bsrc(1)).^2/hvtx(1)^2 - (kit.y(pd(:,i))-Bsrc(2)).^2/hvtx(2)^2 );
end

% solution of background
[u1,g,counts] = mrsolint(Bndl',r,1,nlvl,lvb,bvps,bghs, [],[], FI,SI,sz, Mb,pd);
u1(Bmul(1)+1:Bmul(1):Bndl(1)-Bmul(1),:) = u1(Bmul(1)+1:Bmul(1):Bndl(1)-Bmul(1),:)/2;
u1(:,Bmul(2)+1:Bmul(2):Bndl(2)-Bmul(2)) = u1(:,Bmul(2)+1:Bmul(2):Bndl(2)-Bmul(2))/2;

figure;
surf(kit.x/(hvtx(1)*B0(1)), kit.y/(hvtx(2)*B0(2)) , real(u1) ,'EdgeColor','none');
title('Reference solution');
xlim([0,1]);
ylim([0,1]);
xlabel('x');
ylabel('y');
view(0,90);

%% perturbation
disp('--Solution of perturbation--');
for rid = [4, 16] % root nodes of perturbed subtrees
    rlvl = floor(log2(rid)) + 1;
    nod   = 2^(nlvl-rlvl);
    nodtr = 2^(nlvl-1);
    
    kloc2 = kloc;
    supp = [nsub, [1;1]];
    for il = nod:nod*2-1
        iltr = il + (rid-1)*nod; % index in the global tree
        iatr = iltr - nodtr + 1; % index in the global level
        kloc2{iatr}.val = kloc{iatr}.val/2;
        supp(:,1) = min( supp(:,1), clf(:,iatr) );
        supp(:,2) = max( supp(:,2), clf(:,iatr) );
    end
    
    fprintf('Location: (%d--%d)/%d, (%d--%d)/%d\n', supp(1,1)-1,supp(1,2),nsub(1), supp(2,1)-1,supp(2,2),nsub(2));
    [kit, A2,~,~,~,A2g] = localfem2d(Pele,Pfac,Bvtx,hvtx,bvps,bghs, nlvl,nsub,lvb,clf, kloc2,hvtx);
    
    % new solution update
    [u2, counts] = mrsolupd(u1,rid, nlvl,lvb,bvps,bghs, FI,SI,FE,SE,sz, A1,A2, Mb,pd);
    u2(Bmul(1)+1:Bmul(1):Bndl(1)-Bmul(1),:) = u2(Bmul(1)+1:Bmul(1):Bndl(1)-Bmul(1),:)/2;
    u2(:,Bmul(2)+1:Bmul(2):Bndl(2)-Bmul(2)) = u2(:,Bmul(2)+1:Bmul(2):Bndl(2)-Bmul(2))/2;
    fprintf('  Total time: %f, flops: %.2e\n', sum(counts.time), sum(counts.nflops));
    
    figure;
    surf(kit.x/(hvtx(1)*B0(1)), kit.y/(hvtx(2)*B0(2)) , real(u1+u2) ,'EdgeColor','none');
    ss = sprintf('New solution for updates in (%.2f,%.2f)x(%.2f,%.2f)',...
        (supp(1,1)-1)/nsub(1),supp(1,2)/nsub(1),(supp(2,1)-1)/nsub(2),supp(2,2)/nsub(2));
    title(ss);
    xlim([0,1]);
    ylim([0,1]);
    xlabel('x');
    ylabel('y');
    view(0,90);
end