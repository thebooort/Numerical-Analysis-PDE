
bx = 1;
N=20
N=N+1
m=1000;
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
%vector de coeficientes de los pesoso de los bsplines
a=zeros(1,N+3);

%calculo de los bsplines
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
for n=1:N+3
hold on    
  plot(x,b1(n,1:m+1),'color',rand(1,3))  
hold on
title ("Base de bsplines");
xlabel ("x");
ylabel ("y");
endfor
