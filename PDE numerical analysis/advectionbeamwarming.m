function [x,uint] = advectionbeamwarming(m)
%
% Resuelve u_t + a u_x = 0  en [ax,bx] con condiciones de 
% contorno periodicas empleando el metodo de beam-warmong  modificado 
% con m nodos espaciales interiores.


% Datos del problema
global a
a = -1;           % parametro de la ecuacion
ax = 0;
bx = 1;
tfinal = 1;                % Tiempo maximo
eta= @(x) exp(-600*(x - 0.5).^2); % Cond inicial

% Datos de las discretizaciones
h = (bx-ax)/(m+1);         % h paso espacial
k = 0.4*h;                  % k paso temporal
nu = a*k/h;                % numero de Courant
x = linspace(ax,bx,m+2)';  % Ojo: x(1)=0 y x(m+2)=1
nsteps = round(tfinal / k);    % numero de pasos temp
nplot = 20;       % representamos cada nplot pasos

if -nu/2>1 || -nu/2<0 
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

% En el caso con conciones periodicas se tienen:  
% m+1 incognitas u(2:m+2)  y   u(1) = u(m+2)
I = 2:(m+2);   % indices de las incognitas

% Definicion de condiciones iniciales:
tn = 0;
u0 = eta(x);
u = u0;

% Condiciones de contorno periodicas:
u(1) = u(m+2);   
u(m+3) = u(2);    
u(m+4) = u(3);   % OJO: nodo fantasma para cond periodicas

% representamos condicion inicial:
clf
plot(x,u0)
axis([0 1 -.2 1.2])
title('Condicion inicial a tiempo 0')
input('Pulse <return> para continuar ');

% Bucle temporal:
for n = 1:nsteps
  tnp = tn + k;   % = t_{n+1}

  % Lax-Wendroff vectorial:
  u(I) = u(I) - 0.5*nu*(-3*u(I)+4*u(I+1) - u(I+2)) + 0.5*nu^2 * (u(I+2) + u(I) -2*u(I+1));

  % Condiciones periodicas:
  u(1) = u(m+2);   
  u(m+3) = u(2);   
  u(m+4) = u(3);   % OJO: nodo fantasma para cond periodicas

  % representamos resultados cada nplot pasos:
  if mod(n,nplot)==0 || n==nsteps
    uint = u(1:m+2);  % representamos nodos en el intervalo
    plot(x,uint)
    axis([0 1 -.2 1.2])
    title(sprintf('t = %9.5e  tras %4i pasos con %5i nodos',...
                       tnp,n,m+1))
    if n<nsteps, pause(0.02); end;
  end

  tn = tnp;   % para el siguiente paso temporal
end % del for
uint = u(1:m+2); 
end 
