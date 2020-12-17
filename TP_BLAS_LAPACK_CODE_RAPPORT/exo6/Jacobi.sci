function [x,errf,k]=jacobi(a,b,nit,eps,x0)
   n=size(x0,1);
   x=zeros(n,1);
   k=1;
   errf=zeros(nit,1);
   err=4;
   while k<nit
       x(1,1)=(1/a(1,1))*(b(1)-a(1,2)*x0(2));
   x(n,1)=(1/a(n,n))*(b(n)-a(n,n-1)*x0(n-1));
       for i=2:n-1
           x(i,1)=(1/a(i,i))*(b(i)-a(i,i-1)*x0(i-1)-a(i,i+1)*x0(i+1))
       end
       err=norm(x-x0);
       if err<eps
           break
       end
       errf(k,1)=err;
       k=k+1;
       x0=x
           
   end
endfunction

//Pour tester doit entrer 
//a=[2,1,0;1,2,1;0,1,2]  
//b=[3;4;3]
//x0=[0;0;0]
//[x,errf,k]=jacobi(a,b,40,1e-14,x0)
//tel que nit=nb iterations qu'on veut et la precision qu'on veut=eps=1e-14
