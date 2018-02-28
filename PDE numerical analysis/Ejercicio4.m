%EJERCICIO 4 - BARTOLOM� ORTIZ VISO Y MIGUEL L�PEZ P�REZ
%
%Resolucion de la ecuacion de Laplace sobre el disco:
%
%\Delta u(x)=0 con x \in D={x_1^2+x_2^2\leq 1} 
%
% u(x)=x_1^2  para x\in \mathbb{S}^1=\partial D
%
%transformandola mediante un cambio de variables que nos permite tratar el problema
%en un cuadrado [0,1]\times[0,1] siguiendo las indicaciones del ejercicio
%8 de la relacion 2.
%
%En un paso final se podria adaptar el codigo para compararlo con la solucion 
%analitica exacta,evaluando esta en el mallado y comparando los resultados en valor absoluto
% como en otros ejercicios a entregar.




% DEFINIMOS los datos del problema
a = 0; 
b = 1; 

% DEFINIMOS los datos de la discretizacion
m = 100;     
h = (b-a)/(m+1);

% Puntos del mallado incluidos los contornos
x = linspace(a,b,m+2);  
y = linspace(a,b,m+2); 
[X,Y] = meshgrid(x,y);       
X = X' ;    % De esta forma las coordenadas del nodo (i,j) 
Y = Y' ;    
%
% CONSTRUCCION del SEL: 

e = ones(m,1); %vector de unos
%Generamos matrices C
e1=(e.*((1)./(2*(pi*x(2:m+1)).^2)))(1,:);
C= spdiags([e1'],[0],m,m);
%Creamos lo que llamamos matriz A, mas o menos como en el ejercicio 1, teniendo en cuenta el desfase de  el spdiags
diag=((x(1:m)+(h/2))./(x(1:m)))';
diag2=((x(2:m+1)-(h/2))./(x(2:m+1)))';
diapri=(((x(2:m+1)+(h/2))+(x(2:m+1)-(h/2)))./(x(2:m+1))+((1)./(2*pi*x(2:m+1)).^2));
diapri=diapri';
A= spdiags([diag2 diapri diag],[-1 0 1],m,m);


%donde se va a almacenar la matriz A, la diagonal principal
I = speye(m+1); 
%Ponemos donde se va a almacenar la C, en la posici�n S y T
S = spdiags([e e],[-1 1],m+1,m+1);
T=zeros(m+1,m+1);
T(1,m+1)=1;
T(m+1,1)=1;

%Creamos la matriz A finalmente metiendo las cajas
D = (kron(I,A) + kron(S,C)+kron(T,C)) ./(h^2); 

#T�rminos independientes

#aproximaci�n del origen
aprox=(1/(m+1)).*sum((cos(2.*pi.*x(1:m+1)).^2));
aprox2=ones(1,m).*aprox;

rhs=zeros(m,m+1); #matriz donde guardaremos los t�rminos independientes
a=((x(2:m+1)-h/2)./(x(2:m+1).*(m+1))).*aprox2; 
rhs(:,1)=-a;
b=(cos(2*pi.*x(2:m+1)).^2).*(((x(m+1)+h/2)./x(m+1))); 
rhs(:,m) = -b;
F = reshape(rhs,[],1);

%
% RESOLVEMOS el SEL:
uvec = D\F;  
uvec;
Iint = 2:m+1;       % indices de los nodos interiores en x
Jint = 1:m+1; 
#ponemos el valor de los nodos interiores
usoln(Iint,Jint) = reshape(uvec,m,m+1);
#ponemos los t�rminos de frontera de la soluci�n
usoln(1:m+2,m+2)=cos(2*pi*x(1:m+2)); #condici�n Dirichlet en u(1,theta)
usoln(1:m+2,1)=(1/(m+1))*sum(cos(2.*pi.*x(1:m+1)).^2); #aproximaci�n de la condici�n dirichlet en u(0,theta)
usoln(m+2,1:m+2)=usoln(1,1:m+2); #condici�n peri�dica u(r,0)=u(r,1)
surf(X,Y,usoln)
%pause
%contour3(X,Y,usoln,30)
axis([a b a b])
title('Solucion numerica calculada')