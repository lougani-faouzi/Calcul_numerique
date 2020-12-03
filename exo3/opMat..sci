s=50
rand("seed",s)

n1=25;


A=rand(n1,n1);
B=rand(n1,n1);

tic;
C3=matmat3b(A,B);
t3=toc();

tic;
C2=matmat2b(A,B);
t2=toc();

tic;
C1=matmat1b(A,B);
t1=toc();

tic;
C1k=matmat1b(A,B);
toc

tic;
C=A+B;
toc



n2=50;


A=rand(n2,n2);
B=rand(n2,n2);

tic;
C33=matmat3b(A,B);
tt3=toc();

tic;
C22=matmat2b(A,B);
tt2=toc();

tic;
C11=matmat1b(A,B);
tt1=toc();



n3=75;


A=rand(n3,n3);
B=rand(n3,n3);

tic;
C333=matmat3b(A,B);
ttt3=toc();

tic;
C22=matmat2b(A,B);
ttt2=toc();

tic;
C11=matmat1b(A,B);
ttt1=toc();

xlabel('taille du probleme n')
ylabel('temp exprimé par toc()')
plot2d([n1,n2,n3],[t3,tt3,ttt3],style=[color(255,0,0)])
plot2d([n1,n2,n3],[t2,tt2,ttt2],style=[color(0,255,0)])
plot2d([n1,n2,n3],[t1,tt1,ttt1],style=[color(0,0,255)])




