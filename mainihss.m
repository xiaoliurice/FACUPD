% Main function script for interconnected HSS direct solver for finite
% element systems

addpath('FEM');
addpath('IHSS');

params = struct(...
   'nlvl',  7:2:13,...           % num. of levels of bisection 
   'blvl',  2,...                % num. of bottom levels (no compression)
   'nvtx', [80*[1,2,4,8]; ...    % num. of vertices:             
            80*[1,2,4,8]], ...   %  (rows: x/y; cols: diff. scales)
   'z', 1, ...                   % see below
   'bvps',[2,2,2,2], ...         % boundary condition Up Down Left Right
   'pele', 4,'pfac', 2, ...      % poly order for element and face
   'tol', 1e-6, ...              % HSS rel. compression tolerance
   'low_rank', true,...          % low-rank compression on/off
   'visual', true);              % visualization on/off

 for runId = 1
    % domain partitioning
    nlvl = params.nlvl(runId);
    if(params.low_rank == 1) % select switching level
        slvl = nlvl - params.blvl;
    else
        slvl=0;
    end
    B0   = params.nvtx(:,runId);
    [nsub,lvb,clf, lvh,whb,shb,ihb,lhb,nre] = bisect2(B0,nlvl,slvl);
    [B0, Bvtx] = partition(nsub, 0, B0);
    Bndl = B0*params.pele+1;
    fprintf('--Setup--\nNumber of vertices: global %dx%d, local %dx%d\n',B0,Bvtx);
    fprintf('Number of nodals: %dx%d\n',Bndl);
    if(slvl>0)
        fprintf('Reuse factor: %d/%d=%f\n',nre(1),nre(2), nre(1)/nre(2));
    end
    
    % build local slowness model
    kmax = 2*pi/24; % max. dimensionless wavenumber 2pi/points per wavelength
    kmin = 2*pi/72; % min. wavenumber
    hvtx   = [params.pfac; params.pfac]; % mesh spacing
    w1 = 2; w2 = 4;
    s1 = Bvtx(1)*nsub(1); s2 = Bvtx(2)*nsub(2);
    kloc = cell(1,prod(nsub)); % wavenumber in each loc. subdomain
    for i = 1:size(clf,2)
        % params.z changes the complex argument
        kloc{i} = struct('val',ones(Bvtx'+1)*(kmax*params.z),'d',zeros(size(Bvtx)));
        % border
        if( min([clf(1,i)-1, nsub(1)-clf(1,i), clf(2,i)-1, nsub(2)-clf(2,i) ]) == 0); continue; end
        o2 = Bvtx(2)*(clf(2,i)-1);
        o1 = Bvtx(1)*(clf(1,i)-1);
        for i2 = 1:Bvtx(2)+1
            for i1 = 1:Bvtx(1)+1
                if(i1+o1 < s1*0.45 && (i2+o2-s2*0.5)<s2*0.45 && i1+o1>s1*0.05 && mod(floor((i1+o1)/w2)+floor((i2+o2)/w2),2) == 0)
                    kloc{i}.val(i1,i2) = kmin*params.z;
                end
                if(i1+o1 > s1*0.55 && (i2+o2-s2*0.5)<s2*0.45 && mod(floor((i1+o1)/w1),3) == 0)
                    kloc{i}.val(i1,i2) = kmin*params.z;
                end
            end
        end
    end
    
    % boundary condition
    bghs  = zeros(2,2,3); %1 NtD, 2 ItI (for Helmholtz z=1), 3 RtR (for Poisson z=0)
    bghs(:,:,1) = [0,1;1,0];
    bghs(:,:,2) = [-1i*kmax,1; -1i*kmax,-1];
    bghs(:,:,3) = [kmax,1;kmax,-1];
    
    % generate local finite element systems
    tic;
    [kit, A, M, Mb, pd] = localfem2d(params.pele,params.pfac,Bvtx,hvtx,params.bvps,bghs, nlvl,nsub,lvb,clf, kloc, hvtx);
    t=toc;
    fprintf('Matrix: n %d, nnz %d; FEM time %.1fs\n',prod(Bndl), nnz(A{1})*length(A),t);
    
    % right hand side
    Bsrc = B0.*hvtx.*[0.5;0.5];
    f = zeros(prod(Bndl),1);
    for i = 1: length(f)
        f(i) = exp(- norm([kit.x(i);kit.y(i)]-Bsrc)^2 / hvtx(1)^2 );
    end
    
    disp('--Factorization--');
    [FI,SI,sz, V,rk, counts] = mrfachssint(1,nlvl,slvl,params.tol,lvh,whb,shb,ihb,lhb, lvb,params.bvps,bghs, A,Mb);
    disp(counts);
    fprintf(' time: %f nflops: %.2e storage %.2e\n', sum(counts.time), sum(counts.nflops), sum(counts.storage));
    
    
    disp('--Solution--');
    r = zeros(size(M,1), prod(nsub));
    [u,g,counts] = mrsolhssint(f,r,1,nlvl,slvl,lvb,params.bvps,bghs, [],[],[],[], FI,SI,sz,V,rk, lvh,whb,shb,ihb,lhb, M,Mb,pd);
    disp(counts);
            
    if(params.visual)
        figure;
        title('Wavenumber');
        surf(kit.x,kit.y,real(kit.k),'EdgeColor','none');
        colormap('jet');
        xlim([0,hvtx(1)*B0(1)]);
        ylim([0,hvtx(2)*B0(2)]);
        axis off;
        view(90,90);
        
        figure;
        title('Solution');
        surf(kit.x,kit.y,reshape(real(u),size(kit.x)),'EdgeColor','none');
        colormap('jet');
        xlim([0,hvtx(1)*B0(1)]);
        ylim([0,hvtx(2)*B0(2)]);
        axis off;
        view(90,90);
    end
 end