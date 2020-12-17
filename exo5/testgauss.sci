 A=[3,2;4,5];
 B=[7;6];
 [I,K]=gausskij3b(A,B);
 val=zeros(2);
 val=Usolve(I,K);
 disp(val)


