%% Exame Prático de Recurso Física Computacional — 2017/2018
%% 2 a)
clear all
close all
clc

p = 1.0*10^3;
mi = 1.0*10^-3;
d = 0.010;
G = 800;

t0 = 0; %s
tf = 500; %s
dt = 0.1; %s
t = t0:dt:tf; %s
Nt = length(t);

x0 = 0; 
xf = d; 
dx = 0.0001;
x = x0:dx:xf;
Nx = length(x);

eta = mi*dt/(p*(dx)^2);

u = zeros(Nx, Nt);

A1 = diag(repmat(2/eta + 2,1,Nx-2)); % (2/n + 2) na diagonal

A2 = diag(-ones([1 Nx-3]),1); % sobe um posição relativamente a diagonal
A3 = diag(-ones([1 Nx-3]),-1);% desce uma posição relativamente a diagonal

A = A1 + A2 + A3;
b = zeros(Nx-2,1);

const = 2/eta - 2;
c2 = dt/(eta*p);


for j = 1: Nt -1
    b = u(1:Nx-2,j) + const * u(2:Nx-1,j) + u(3:Nx,j)+ c2 *(G+G); 
    b(1) = b(1) + u(1,j+1);     % o que esta a vermelho no slide 10 - Aula 6 (2020/21)
    b(Nx-2) = b(Nx-2) + u(Nx,j+1);  % o que esta a vermelho no slide 10 - Aula 6 (2020/21)
    u(2:Nx-1,j+1) = linsolve(A,b);

end

figure(1)
mesh(t,x,u), xlabel('t(s)'), ylabel('x(cm)'), zlabel('T(ºC)')
figure(2)
contourf(x,t,u'), xlabel('x(cm)'), ylabel('t(s)')
h = colorbar; % barra de cores
set(get(h,'label'),'string','T(ºC)'); % label da barra de cores

%% 1 b)
clear all
close all
clc

p = 1.0*10^3;
mi = 1.0*10^-3;
d = 0.010;

t0 = 0; %s
tf = 500; %s
dt = 0.1; %s
t = t0:dt:tf; %s
Nt = length(t);

x0 = 0; 
xf = d; 
dx = 0.0001;
x = x0:dx:xf;
Nx = length(x);

G = 800*tanh(t/10);

eta = mi*dt/(p*(dx)^2);

u = zeros(Nx, Nt);

A1 = diag(repmat(2/eta + 2,1,Nx-2)); % (2/n + 2) na diagonal

A2 = diag(-ones([1 Nx-3]),1); % sobe um posição relativamente a diagonal
A3 = diag(-ones([1 Nx-3]),-1);% desce uma posição relativamente a diagonal

A = A1 + A2 + A3;
b = zeros(Nx-2,1);

const = 2/eta - 2;
c2 = dt/(eta*p);

for j = 1: Nt -1
    b = u(1:Nx-2,j) + const * u(2:Nx-1,j) + u(3:Nx,j)+ c2 *(G(j)+G(j+1)); 
    b(1) = b(1) + u(1,j+1); 
    b(Nx-2) = b(Nx-2) + u(Nx,j+1);
    u(2:Nx-1,j+1) = linsolve(A,b);

end

figure(1)
mesh(t,x,u), xlabel('t(s)'), ylabel('x(cm)'), zlabel('T(ºC)')
figure(2)
contourf(x,t,u'), xlabel('x(cm)'), ylabel('t(s)')
h = colorbar; % barra de cores
set(get(h,'label'),'string','T(ºC)'); % label da barra de cores

%% 1 c)
clear all
close all
clc

p = 1.0*10^3;
mi = 1.0*10^-3;
d = 0.010;

t0 = 0; %s
tf = 500; %s
dt = 0.1; %s
t = t0:dt:tf; %s
Nt = length(t);

x0 = 0; 
xf = d; 
dx = 0.0001;
x = x0:dx:xf;
Nx = length(x);

G = 800*tanh(t/10);

eta = mi*dt/(p*(dx)^2);

u = zeros(Nx, Nt);
% Condição inicial
u(:,1) = 5*x/d;

A1 = diag(repmat(2/eta + 2,1,Nx-2)); % (2/n + 2) na diagonal

A2 = diag(-ones([1 Nx-3]),1); % sobe um posição relativamente a diagonal
A3 = diag(-ones([1 Nx-3]),-1);% desce uma posição relativamente a diagonal

A = A1 + A2 + A3;
b = zeros(Nx-2,1);

const = 2/eta - 2;
c2 = dt/(eta*p);

% for j = 1: Nt -1
%     b = u(1:Nx-2,j) + const * u(2:Nx-1,j) + u(3:Nx,j)+ c2 *(G(j)+G(j+1)); 
%     b(1) = b(1) + u(1,j+1); 
%     b(Nx-2) = b(Nx-2) + u(Nx,j+1);
%     u(2:Nx-1,j+1) = linsolve(A,b);
% 
% end

for n = 1:Nt-1
    for i = 1: Nx-2
        b(i) = u(i,n) + const*u(i+1,n) + u(i+2,n)+ c2 *(G(n)+G(n+1));
    end
    b(1) = b(1) + u(1,n+1); 
    b(Nx-2) = b(Nx-2) + u(Nx,n+1);
    u(2:Nx-1,n+1) = linsolve(A,b);
end

figure(1)
mesh(t,x,u), xlabel('t(s)'), ylabel('x(cm)'), zlabel('T(ºC)')
figure(2)
contourf(x,t,u'), xlabel('x(cm)'), ylabel('t(s)')
h = colorbar; % barra de cores
set(get(h,'label'),'string','T(ºC)'); % label da barra de cores