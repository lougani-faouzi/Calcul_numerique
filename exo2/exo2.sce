function exo2(n,ex)
     format("e",ex);
     A = rand(n,n);
     disp("A=",A);
     xex = rand(n,1);
     disp("xex=",xex);
     b = A*xex;
     disp("b=",b);
     x=A\b;
     disp("x=",x);
     frelres = norm(x-xex)/norm(xex); //erreur avant
     disp("frelres=",frelres);
     brelres = norm(b-A*x)/norm(b); //erreur arrière
     disp("brelres=",brelres);
     capa=cond(A);
     disp("cap=",capa);
     borne=cond(A)*brelres;
     disp("borne=",borne);
endfunction


    
    
