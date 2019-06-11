function x = mrsolvec(b, sizeu, Bmul, nlvl,lvb,bvps,bghs, FI,SI,sz,Mb,pd)
    r = zeros(size(pd));
    for i = 1:size(pd,2)
        r(:,i) = b(pd(:,i));
    end
    x = mrsolint(sizeu,r,1,nlvl,lvb,bvps,bghs, [],[],FI,SI,sz, Mb,pd);
    x(Bmul(1)+1:Bmul(1):sizeu(1)-Bmul(1),:) = x(Bmul(1)+1:Bmul(1):sizeu(1)-Bmul(1),:)/2;
    x(:,Bmul(2)+1:Bmul(2):sizeu(2)-Bmul(2)) = x(:,Bmul(2)+1:Bmul(2):sizeu(2)-Bmul(2))/2;
    x = reshape( x, size(b) );
end