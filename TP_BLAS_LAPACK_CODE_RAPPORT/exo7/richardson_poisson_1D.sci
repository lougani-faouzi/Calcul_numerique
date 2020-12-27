function[x,relres,resvec,it]=richardson_poisson_1d(A,b,tol,maxit,x0,alpha)
    it=0;
    res=b-A*x0;
    relres=norm(b-A*x0)/norm(b);
    while (relres>tol)&(it<maxit)
    it=it+1;
    x=x0+alpha;
    res=b-A*x0;
    relres=norm(b-A*x0)/norm(b);
    resvec(it)=relres;
    x0=x;
    end   
endfunction
