function [B1, Bi] = partition(nsub, novl, B0)
% B0:   initial grid
% nsub: number of blocks on each direction
% novl: number of overlapping grids between adjacent subdomains
% B1:   modified grid so that each subdomain is equal in size
% Bi:   size of each subdomain

nd  = length(nsub);

B1 = B0;
Bi = B0;

for i = 1:nd
    l = (nsub(i)-1)*novl + B0(i);
    k = mod(l,nsub(i));
    if (k~=0)
        B1(i) = B0(i) + (nsub(i)-k);
        l     = l     + (nsub(i)-k);
    end
    Bi(i) = l/nsub(i);
end

end
