ax = 0;
bx = 1;
N=90;
N=N+1;
m=100;
k = (bx-ax)/(N); 
h=k;   
nodos = linspace(ax,bx,N)';     
x = linspace(ax,bx,m+1)';
m=100;
t=1;
f = @(x) sin(pi*x);
final = @(x,y) exp((-pi^2)*y).*sin(pi*x);
%vector de coeficientes de los pesos de los b-splines
a=zeros(1,N+3);

#cálculo de la base de b-splines
b1=zeros(N+3,m+1);
for n = 1:3
    if (n==1)
    for j=1:m+1
    b1(1,j)=B_11(x(j)/k);
    b1(N+3,m+2-j)=b1(1,j);
    endfor
    
    elseif (n==2)
    for j=1:m+1
    b1(2,j)=B_0(x(j)/k);
    b1(N+2,m+2-j)=b1(2,j);
    endfor
    
    
    elseif (n==3)
    for j=1:m+1
    b1(3,j)=B_1(x(j)/k);
    b1(N+1,m+2-j)=b1(3,j);
    endfor
    endif
    endfor
for n=4:N
   for j=1:m+1
   b1(n,j)=B_n((n-4),x(j)/k) ;  
   endfor
   endfor 

%for n=1:N+3
%hold on    
%  plot(x,b1(n,1:m+1),'color',rand(1,3))  
%hold on
%title ("Base de bsplines");
%xlabel ("x");
%ylabel ("y");
%endfor




A=zeros(N-1,N-1);
B=zeros(N-1,N-1);

%Diagonal principal de A
diagop1=zeros(1,N-1);
for j=1:N-1
  diagop1(j)=2/3
  endfor
%Subdiagonales de A
subdiago1=zeros(1,N-1);
for j=1:N-1
  subdiago1(j)=1/6
  endfor

%Diagonal principal de B
diagop2=zeros(1,N-1);
for j=1:N-1
  diagop2(j)=-2/k^2;
  endfor
  
%Subdiagonales de B
subdiago2=zeros(1,N-1);
for j=1:N-1
  subdiago2(j)=1/k^2;
  endfor
  
%Reorientamos los vectores diagonal
diagop1=diagop1'
diagop2=diagop2'
subdiago1=subdiago1'
subdiago2=subdiago2'
%creamos las matrices A y B
A = spdiags([subdiago1 diagop1 subdiago1],[-1 0 1],N-1,N-1);
B = spdiags([subdiago2 diagop2 subdiago2],[-1 0 1],N-1,N-1);
%Primer vector a_0
a_0=zeros(1,N-1)';

%vector de evaluaciones de g para calcular el a_0
indepe=zeros(1,N-1);
for i=1:N-1
    indepe(i)=f((i-1)*k);
endfor
indepe=indepe';
a_0=A\indepe;
for n=3:N+1
    a(n)=a_0(n-2)
    endfor

a(1)=0;
a(N+3)=0;
a(2)=(1/3)*(a(3));
a(N+2)=(1/3)*(a(N+1));

for n=1:N+3   
  b1(n,1:m+1)=a(n)*b1(n,1:m+1) ;
endfor
solucion=zeros(1,m+1);

for n=1:m+1
aux=0; 
     for j=1:N+3
         aux=aux+b1(j,n) ; 
         endfor
  sol(n)=aux ;
endfor
plot(x,sol,'r')
hold on
%evolucion temporal:
tini=0;
tfin=0.1;
pasos=10000;
tep=(0.1/pasos);
t=1;

for t=1:7000
    a_0=a_0.+tep*((inv(A)*B)*a_0);
    t=1+t;
    endfor

a(1)=0;
a(N+3)=0;
a(2)=(1/3)*(a(3));
a(N+2)=(1/3)*(a(N+1));

for n=1:N+3   
  b1(n,1:m+1)=a(n)*b1(n,1:m+1) ;
endfor
solucion=zeros(1,m+1);

for n=1:m+1
aux=0; 
     for j=1:N+3
         aux=aux+b1(j,n) ; 
         endfor
  sol(n)=aux ;
endfor
plot(x,sol,'b')
hold on
legend ("aproximacion inicial","aproximacion en 0.01");