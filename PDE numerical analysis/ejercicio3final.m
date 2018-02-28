% Ejercio 3 relacion PVF's
%
% Resuelve el problema de valores en la frontera
% \Delta u=2\pi^2sin(\pi y)cos(\pi x)
%
%con (x,y)\in [0,1]\times[0,1]  
% Y valores iniciales
%u(0,y)=-u(1,y)=sin(\pi y),u(x,0)=-u(x,1)=0
%
% Se emplea un esquema de 9 puntos: discretización del operador la-
%placiano en dimensión dos con una molecula computacional de 9 puntos
%
%El programa resuelve directamente el PVF y representa la solucion

% DEFINIMOS los datos del problema
a = 0; 
b = 1; 
% funcion f(x,y) (evaluable sobre matrices)
f = @(x,y) 2*pi^2*cos(pi*x).*sin(pi*y);      
%
% DEFINIMOS los datos de la discretizacion
m = 10;      
h = (b-a)/(m+1);
% Puntos del mallado incluidos los contornos
x = linspace(a,b,m+2);  
y = linspace(a,b,m+2); 
[X,Y] = meshgrid(x,y);       
X = X' ;     
Y = Y' ;    
%
% CONSTRUCCION del SEL: 
I = speye(m);  
e = ones(m,1); 
T = spdiags([4*e -20*e 4*e],[-1 0 1],m,m);  
S = spdiags([e e],[-1 1],m,m);
U = spdiags([e 4*e e],[-1 0 1],m,m); 
A = (kron(I,T) + kron(S,U)) /(6* h^2); 
% definimos el termino independiente que tiene en cuenta 
% los valores en la frontera.
Iint = 2:m+1 ;      
Jint = 2:m+1 ;
Iext1 = 1:m;
Iext2 = 3:m+2;
Jext1 = 1:m;
Jext2 = 3:m+2;      
Xint = X(Iint,Jint);       
Yint = Y(Iint,Jint);
rhs = f(Xint,Yint);        
% funcion que determina el valor en la frontera.
% Aunque se evalua en todos los nodos, solo se usaran 
% los de la frontera a continuacion. De nuevo tomamos el valor de la funcion 
%solucion para m´as exactitud, en caso de tener que usar los datos de las condiciones
%iniciales basta con evaluar las funciones alli definidas en los terminos especificos
usoln = -cos(pi*X).*sin(pi*Y);             
% se incorpora en rhs los valores en la frontera calculados con la formula
rhs(:,1) = rhs(:,1) - 4*usoln(Iint,1)/(6*h^2)-usoln(Iext1,1)/(6*h^2)-usoln(Iext2,1)/(6*h^2);
rhs(:,m) = rhs(:,m) - 4*usoln(Iint,m+2)/(6*h^2)-usoln(Iext1,m+2)/(6*h^2)-usoln(Iext2,m+2)/(6*h^2);
rhs(1,:) = rhs(1,:) - 4*usoln(1,Jint)/(6*h^2)-usoln(1,Jext2)/(6*h^2)-usoln(1,Jext1)/(6*h^2);
rhs(m,:) = rhs(m,:) - 4*usoln(m+2,Jint)/(6*h^2)-usoln(m+2,Jext2)/(6*h^2)-usoln(m+2,Jext1)/(6*h^2);
%En la realizacion de rhs veo que los terminos extremos estan repetidos dos veces
%divido entre dos para solucionar el problema
rhs(1,m) = rhs(1,m)/2;
rhs(m,m) = rhs(m,m)/2;
rhs(1,1) = rhs(1,1)/2;
rhs(m,1) = rhs(m,1)/2;

% Convertimos la matriz rhs en un vector columna  
F = reshape(rhs,m*m,1);  
%
% RESOLVEMOS el SEL:
uvec = A\F;  
%
% REPRESENTAMOS la solucion numerica
usoln(Iint,Jint) = reshape(uvec,m,m);
surf(X,Y,usoln);
axis([a b a b]);
title('Solucion numerica calculada');