function [t,y]=eulerimplicito(A,t0,tfin,u0,F,m)
%Resuelve un sistema implícito dado
%
% Devuelve: 
% t vector con los nodos
% y valores aproximados: y(:,i) aproximacion  de u(t(i)).

%A partir de el primer gap ya no se muestra en el help 
% Inicializamos variables
k= (tfin-t0)/(m);          
t=linspace(t0,tfin,m+1);      
y(:,1)=u0;
I=eye(m);
M=I-k*A;
for n=1:m
     y(:,n+1)=M\(y(:,n)+k*(F));
end;
y=[F,y,F];


