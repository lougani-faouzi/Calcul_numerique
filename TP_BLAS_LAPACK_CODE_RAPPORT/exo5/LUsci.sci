function [l,d] = facto(a,b,c)
    n = length(b);
    d(1) = b(1);
    for i=2:n
        l(i-1) = a(i-1)/d(i-1);
        d(i)   = b(i)-l(i-1)*c(i-1);
    end
endfunction
