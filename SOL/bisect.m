function [nsub, lvb, clf] = bisect(B0, nlvl)
% nsub:  number of blocks on each direction
% lvb:   level of boundary
% clf:   coordinate of leaf
nsub   = ones(size(B0));
lvb    = cell(nlvl,1);
lvb{1} = ones(4,1);
nn = 1;
clf = [1;1];
for lvl = 2: nlvl
    [~,id] = max(B0./nsub);
    nsub(id) = nsub(id)*2;
    lvb{lvl} = ones(4,nn*2);
    lvb{lvl}(:     ,1:2:nn*2-1) = lvb{lvl-1};
    lvb{lvl}(:     ,2:2:nn*2)   = lvb{lvl-1};
    lvb{lvl}(id*2  ,1:2:nn*2-1) = lvl;
    lvb{lvl}(id*2-1,2:2:nn*2)   = lvl;
    clf1(:, 1:2:nn*2-1) = clf;
    clf1(:, 2:2:nn*2)   = clf;
    clf1(id,1:2:nn*2-1) = clf1(id,1:2:nn*2-1)*2-1;
    clf1(id,2:2:nn*2)   = clf1(id,1:2:nn*2-1)+1;
    nn = nn*2;
    clf = clf1;
end
lvb = horzcat(lvb{:});
end