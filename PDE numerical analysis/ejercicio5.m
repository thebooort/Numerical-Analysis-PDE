function [x,uint] = ejercicio5(m)
%
% Resuelve u_t + a u_x = f  en [ax,bx] con condiciones de 
% contorno dirichlet empleando el metodo de Lax-Wendroff modificado
% con m nodos espaciales interiores.

% Datos del problema
global a
a = -4;           % parametro de la ecuacion
ax = -10;
bx = 2;
tfinal = 1;                % Tiempo maximo
eta= @(x) 1-exp(-1*(x).^2); % Cond inicial
f=@(t,x) exp(t).-exp(t.-(x.+4*t).^2);
% Datos de las discretizaciones
h = (bx-ax)/(m+1);         % h paso espacial
k = 0.09*h;                  % k paso temporal
nu = a*k/h;                % numero de Courant
x = linspace(ax,bx,m+2)';  % Ojo: x(1)=0 y x(m+2)=1
nsteps = round(tfinal / k);    % numero de pasos temp
nplot = 20;       % representamos cada nplot pasos

if abs(nu)>1 
  disp(' ')
  disp(sprintf('OJO! No se verifica la cond CFL'))
  disp(' ')
end

% Comprobamos si el ultimo paso llega a tfinal
if abs(k*nsteps - tfinal) > 1e-5
  disp(' ')
  disp(sprintf('OJO *** k no divide a tfinal, k = %9.5e',k))
  disp(' ')
end

% En el caso con conciones dirichlet tenemos
I = 2:(m+1);   % indices de las incognitas
% Definicion de condiciones iniciales:
tn = 0;
u0 = eta(x);
u = u0;

% Condiciones de contorno dirichlet:
u(1) = e^tn;   
u(m+2) = e^tn;   
% representamos condicion inicial:
clf
plot(x,u0)
axis([-10 2 -2 6])
title('Condicion inicial a tiempo 0')
input('Pulse <return> para continuar ');
% Bucle temporal:
for n = 1:nsteps
  tnp = tn + k;   % = t_{n+1}

  % Lax-Wendroff vectorial:
  u(I) = u(I) - 0.5*nu*(u(I+1) - u(I-1)) + 0.5*nu^2 * (u(I-1) - 2*u(I) + u(I+1)) + (k*(exp(tn)-exp(tn-(x(I).+4*tn).^2))) + (((-a*k^2)/(4*h))*((exp(tn)-exp(tn-(x(I+1)+4*tn).^2)).-((exp(tn)-exp(tn-(x(I-1)+4*tn).^2)))))+((k/2)*(((exp(tn+k)-exp(tn+k-(x(I).+(4*tn+4*k)).^2))).-((exp(tn)-exp(tn-(x(I)+4*tn).^2)))));

  % Condiciones contorno:
  u(1) = e^tnp;   
  u(m+2) = e^tnp;   

  % representamos resultados cada nplot pasos:
  if mod(n,nplot)==0 || n==nsteps
    uint = u(1:m+2);  % representamos nodos en el intervalo
    plot(x,uint)
    axis([-10 2 -2 6])
    title(sprintf('t = %9.5e  tras %4i pasos con %5i nodos',...
                       tnp,n,m+1))
    if n<nsteps, pause(0.2); end;
  end

  tn = tnp;   % para el siguiente paso temporal
end % del for
uint = u(1:m+2); 
end 
