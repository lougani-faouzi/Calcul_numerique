function A= mylu3b(A)
    n=size(A,1);
    
    for k=1:n-1
        for i=k+1:n
            A(i,k)=A(i,k)/A(k,k);
        end
          for i=k+1:n
            for j=k+1:n
                A(i,j)=A(i,j)-A(i,k)*A(k,j);
            end 
        end
        
    end
endfunction




function A= mylu3b1(A)
   
    n=size(A,1);
    
    for k=1:n-1
        
            A(k+1:n,k)=A(k+1:n,k)/A(k,k);
            A(k+1:n,k+1:n)=A(k+1:n,k+1:n)-A(k+1:n,k)*A(k,k+1:n);
     end 

endfunction



function A= mylu(A)
   
    n=size(A,1);
    // on ajoute la partie pivot 
    for j=1:n-1
        for i=1:n-1
    
         A(i,:)=A(i,:)-A(i,j)/A(j,j)*A(j;:)
       
        end 
    end
    // fin parite pivot 
    for k=1:n-1
        
            A(k+1:n,k)=A(k+1:n,k)/A(k,k);
            A(k+1:n,k+1:n)=A(k+1:n,k+1:n)-A(k+1:n,k)*A(k,k+1:n);
    end 

endfunction


