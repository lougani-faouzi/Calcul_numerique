function [A,b]=gausskij3b(A,b)
    n=size(A,1);
    m=eye(n,n);
    
    for k=1:n-1
        for i=k+1:n
            m(i,k)=A(i,k)/A(k,k);
            b(i)=b(i)-m(i,k)*b(k);
            for j=k+1:n
                A(i,j)=A(i,j)-m(i,k)*A(k,j);
            end 
        end
        
    end
endfunction



function x=Usolve(U,b)
    n=size(U,1);
    x=zeros(n,1);
    x(n)=b(n)/U(n,n);
    for i=n-1:-1:1
       
        x(i)=(b(i)-U(i,(i+1):n)*x(i+1:n))/U(i,i);
    end
endfunction
