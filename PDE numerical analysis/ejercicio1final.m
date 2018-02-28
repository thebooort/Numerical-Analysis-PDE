%Ejercicio 1 de la 1º relacion a entregar:
%Estudio de la ecuacion: 
%d/dx(e^xu'(x))-u(x)=x^2, x\in(-1,1), u(-1)=2,u(1)=3
%A traves de la formula de aproximacion del ejercicio 7:
% aproximacion del termino con derivadas mediante diferencias finitas centradas y paso h/2.

%Definimos los datos del problema
a = -1;     
b = 1; 
alpha = 2;
beta = 3;
% definimos las funciones a utilizar
f = @(x) x.^2; 
g= @(x) exp(x);               
% DEFINIMOS los datos de la discretizacion
m = 100;                   
h = (b-a)/(m+1); 

% Puntos del mallado incluidos los contornos
x = linspace(a,b,m+2);   
%
% CONSTRUCCION del SEL
diag=exp(x(1:m)+(h/2))';
diag2=exp(x(2:m+1)+(h/2))';
diapri=(-(g(x(2:m+1)+(h/2))+g(x(2:m+1)-(h/2)))-h^2);
diapri=diapri';
A= spdiags([diag2 diapri diag],[-1 0 1],m,m);
A= A/h^2;
% definimos el termino independiente que tiene en cuenta 
% los valores en la frontera
xint=(x(2:m+1))';%el apostrofe significa transponer
rhs=f(xint);%rhs = right hand side/termino de la derecha
rhs(1) = rhs(1) -alpha*exp((x(2)-h/2))/h^2;
rhs(m) = rhs(m)- beta*exp((x(m+1)+h/2))/h^2;
%Resolucion del SEL
uvec = A\rhs;
% REPRESENTAMOS la solucion numerica
usoln = [alpha; uvec; beta];%añadimos las condiciones de contorno  y los nodos intermedios
clf;
plot(x,usoln) %generamos puntos aislados
