%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% resolucion de la ecuacion del calor con evolucion temporal.
% Aplica el metodo de lineas al problema de evolucion
% u_t = u_{xx}  , x\in[a,b] t \in [0,T],  u(a,t)=0, u (b,t)=0 
% Discretiza en espacio mediante diferencias finitas centradas
% Ojo!!! Llama a la funci�n Eulerimplicito
% Datos de la ecuacion e intervalo espacial
clear('A');
global A 
a = 0; b = 1;  
alpha = 0; beta =0; 
% datos de la discretizacin
m = 50;                   % numero de nodos interiores
h = (b-a)/(m+1);
%
% Defino la matriz 
e = ones(m,1);
A = spdiags([e -2*e e],[-1 0 1],m,m);			
A= A/h^2;
%condicion inicial
F=zeros(m,1);
% dato inicial para la resolucion          
x=linspace(a,b,m+2);
xint=x(2:m+1);
u=sin(pi*xint)';
plot(xint,u)    
% Resolvemos mediante Euler implicito
[t,sol]=eulerimplicito(A,0,2,u,F,m);
% Representamos la evolucion de la soluciones numerica y exacta  
solexact =@(x,t) exp(-pi^2*t)*sin(pi*x);
clf
hold off
for i=0:1:m
plot(x,[0;sol(:,i+1);0],x,solexact(x,t(i+1)),pause(0.4));
endfor