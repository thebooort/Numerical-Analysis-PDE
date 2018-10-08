ax = 0;
bx = 1;
N=1000
N=N+1
m=100;
k = (bx-ax)/(N); 
h=k;        
x = linspace(ax,bx,m+1)';
m=100;
p = @(x) 1.+x.*0;
q = @(x) 0.+x.*0;
f = @(x) (-e.^(x-1.))-1.
final= @(x) x.*(1.-e.^(x-1))
%condiciones iniciales
alpha_1=0;
beta_1=final(1.1);
%vector de coeficientes de los pesos de los b-splines
a=zeros(1,N+3);

%cálculo de los b-splines
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
%print de la base de b-splines
%for n=1:N+3
%hold on    
 % plot(x,b1(n,1:m+1),'color',rand(1,3))  
%hold on
%title ("Base de bsplines");
%xlabel ("x");
%ylabel ("y");
%endfor

%calculo de la matriz del sistema
matriz_1=zeros(N+1,N+1);

diago=zeros(1,N+1);
for j=3:N-1
  diago(j)=-2/h^2.+(2/3.)*q(j);
  endfor
diago(1)=-9./h^2+3.*p(1)/h
diago(2)=-5./(2*h^2)+(1./4)*(p(2)/h)+(7./12)*q(2)
diago(N)=-5./(2*h^2)-(1./4)*(p(N)/h)+(7./12)*q(N)
diago(N+1)=-9./h^2-3.*p(N+1)/h

diago1=zeros(1,N+1);
for j=4:N
  diago1(j)=1/h^2+1./2*h+(1/6)*q(j);
  endfor
diago1(2)=3./h
diago1(3)=1./(h^2)+(1./2)*(p(2)/h)+(1./6)*q(2)
diago1(N+1)=3./(2*h^2)+(3./4)*(p(N)/h)+(1./4)*q(N)


diago2=zeros(1,N+1);

for j=2:N-1
  diago2(j)=1/h^2-1./2*h+(1/6)*q(j);
  endfor
diago2(1)=3./(2*h^2)-(3./4)*(p(2)/h)+(1./4)*q(2)
diago2(N)=3./h^2

diago=diago'
diago2=diago2'
diago1=diago1'

matriz_1 = spdiags([diago2 diago diago1],[-1 0 1],N+1,N+1);
%calculo de los terminos independientes
indepe=zeros(1,N+1)';
for i=2:N
    indepe(i)=f((i-1)*k);
endfor
indepe(1)=f(0)-alpha_1*(6/h^2-3*p(0)/h+q(0));
indepe(N+1)=f(N*k)-alpha_1*(6/h^2-3*p(N)/h+q(N));
%Resolucion del sistema y valores extra
a=matriz_1\indepe;

a(1)=alpha_1;
a(N+3)=beta_1;

%valores finales de la aporximacion
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
plot(x,final(x))
hold on 
plot(x,sol,'r')
title ("Funcion final y aproximacion por b-splines");
xlabel ("x");
ylabel ("y");
legend ("Funcion exacta","aproximacion");