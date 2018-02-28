% Ejercio 2 relacion PVF's
%
% Resuelve el problema de valores en la frontera
% \Delta u=2\pi^2sin(\pi x)cos(\pi y)
%
%con (x,y)\in [0,1]\times[0,1]  
% Y valores iniciales
%u(x,0)=-u(x,1)=sin(\pi x),u(0,y)=-u(1,y)=0
%
% Se emplea un esquema en diferencias en 5 puntos.
% Modificada de  
%http://www.amath.washington.edu/~rjl/fdmbook/chapter3  (2007)
%
%El programa resuelve directamente el PVF y estudia los errores cometidos
%a traves de la diferencia con la solucion analitica. Como resultado
%muestra una representacion del valor absoluto de los errores cometidos

% DEFINIMOS los datos del problema
a = 0; 
b = 1; 
% funcion f(x,y) 
f = @(x,y) 2*pi^2*sin(pi*x).*cos(pi*y);      
%
% DEFINIMOS los datos de la discretizacion
m = 100;     
h = (b-a)/(m+1);
% Puntos del mallado incluidos los contornos
x = linspace(a,b,m+2);  
y = linspace(a,b,m+2); 
[X,Y] = meshgrid(x,y);       
X = X';    
Y = Y';     

% CONSTRUCCION del SEL: 
I = speye(m); 
e = ones(m,1); 
T = spdiags([e -4*e e],[-1 0 1],m,m);  
S = spdiags([e e],[-1 1],m,m); 
A = (kron(I,T) + kron(S,I)) / h^2;

% definimos el termino independiente que tiene en cuenta 
% los valores en la frontera.
Iint = 2:m+1;       
Jint = 2:m+1;       
Xint = X(Iint,Jint);       
Yint = Y(Iint,Jint);
rhs = f(Xint,Yint);  
       
% funcion que determina el valor en la frontera.
% Aunque se evalua en todos los nodos, solo se usaran 
% los de la frontera a continuacion. Creo que es mas eficiente y exacto de esta forma
%Sin embargo dado que tenemos las condiciones iniciales, tambien podrian 
%utilizarse estas para calcular estos valores
usoln = -cos(pi*Y).*sin(pi*X);               
% se incorpora en rhs los valores en la frontera
rhs(:,1) = rhs(:,1) - usoln(Iint,1)/h^2;
rhs(:,m) = rhs(:,m) - usoln(Iint,m+2)/h^2;
rhs(1,:) = rhs(1,:) - usoln(1,Jint)/h^2;
rhs(m,:) = rhs(m,:) - usoln(m+2,Jint)/h^2;
% Convertimos la matriz rhs en un vector columna  
F = reshape(rhs,m*m,1); 

% RESOLVEMOS el SEL:
uvec = A\F;  
%
% REPRESENTAMOS los errores cometidos mediante la resta de los valores calculados
%y los valores exactos
usoln(Iint,Jint) = reshape(uvec,m,m);
usoln=usoln.-(-cos(pi*Y).*sin(pi*X));
A=abs(usoln);
surf(X,Y,A)
axis([a b a b])
title('Errores (en valors absoluto) cometidos')


