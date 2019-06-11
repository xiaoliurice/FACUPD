function b = localmatvec(x, A, pd)
    b = zeros(size(x));
    for i = 1:length(A)
       b(pd(:,i)) = b(pd(:,i)) + A{i}*x(pd(:,i));
    end
end