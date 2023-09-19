%% 2o Teste Prático de Avaliação Discreta Física Computacional — 2017/2018
%% 1 a)
clear all
close all
clc

% Constantes
D = 0.02;
m = 5;
L = 1;

x0 = 0;
xf = L;
t0 = 0;
tf = 5;

dx = 0.01;
dt = 0.001;

x = x0:dx:xf;
t = t0:dt:tf;

Nx = length(x);
Nt = length(t);

u = zeros(Nx,Nt);

% Condições iniciais
u(:,1) = 0.5;
u((x<0.6),1) = 0;
u((x>0.8),1) = 0;

% Condições fronteira
u(1,:) = 0;
u(Nx,:) = 0;

% Método de Euler
for j = 1:Nt-1
    for i = 2:Nx-1
        u(i,j+1) = u(i,j) + D*dt/(dx^2) * (- 2*u(i,j)+ u(i-1,j) + u(i+1,j)) + m * dt * (u(i,j)*(1-u(i,j)));
    end
end

figure(1)
mesh(t,x,u), xlabel('t'), ylabel('x'), zlabel('u')
figure(2)
contourf(x,t,u'), xlabel('x'), ylabel('t'), zlabel('u')
h = colorbar; % barra de cores
set(get(h,'label'),'string','Evolução do gene'); % label da barra de cores

%% 1 b)
clear all
close all
clc

% Constantes
D = 0.02;
m = 5;
L = 1;

x0 = 0;
xf = L;
t0 = 0;
tf = 3;

dx = 0.01;
dt = 0.001;

x = x0:dx:xf;
t = t0:dt:tf;

Nx = length(x);
Nt = length(t);

u = zeros(Nx,Nt);

% Condições iniciais
u(1:end-1,1) = 0.0;

% Condições fronteira
u(1,:) = 0;
u(Nx,:) = 1;

% Método de Euler
for j = 1:Nt-1
    for i = 2:Nx-1
        u(i,j+1) = u(i,j) + D*dt/(dx^2) * (- 2*u(i,j)+ u(i-1,j) + u(i+1,j)) + m * dt * (u(i,j)*(1-u(i,j)));
    end
end

figure(1)
mesh(t,x,u), xlabel('t'), ylabel('x'), zlabel('u')
figure(2)
contourf(x,t,u'), xlabel('x'), ylabel('t'), zlabel('u')
h = colorbar; % barra de cores
set(get(h,'label'),'string','Evolução do gene'); % label da barra de cores

% Olhando para o contourplot vemos que existem bandas
% paralelas da mesma cor, o que significa que entre esses pontos, estas
% bandas paralelas viajam a velocidades constantes


[row, col] = find(u > 0.1-0.001  & u <0.1 + 0.001);
figure(3)
plot(t(col(50:200)),x(row(50:200)))
xlabel('t')
ylabel('x')
p = polyfit(t(col(50:200)),x(row(50:200)),1);

vel = abs(p(1));
v_analitica = 2*sqrt(D*m);

erro = (abs(v_analitica-vel)/v_analitica)*100;
disp(['Velocidade obtida: ',num2str(vel)])
disp(['Velocidade analitica: ',num2str(v_analitica)])
disp(['Erro: ',num2str(erro),' %'])

%% 1 c)
clear all
close all
clc

% Constantes
D = 0.02;
m = 5;
L = 1;

x0 = 0;
xf = L;
t0 = 0;
tf = 3;

dx = 0.01;
dt = 0.001;

x = x0:dx:xf;
t = t0:dt:tf;

Nx = length(x);
Nt = length(t);

u = zeros(Nx,Nt);

% Condições iniciais
u(:,1) = 0;

% Condições fronteira
u(1,:) = 0;
u(Nx,:) = 1;

% Método de Euler
for j = 1:Nt-1
    u(1,j+1) = u(1,j) + (D*dt/dx^2)*(-2*u(1,j)+2*u(2,j)) + m*dt*(u(1,j)*(1-u(1,j)));
    for i = 2:Nx-1
        u(i,j+1) = u(i,j) + D*dt/(dx^2) * (- 2*u(i,j)+ u(i-1,j) + u(i+1,j)) + m * dt * (u(i,j)*(1-u(i,j)));
    end
end

figure(1)
mesh(t,x,u), xlabel('t'), ylabel('x'), zlabel('u')
figure(2)
contourf(x,t,u'), xlabel('x'), ylabel('t'), zlabel('u')
h = colorbar; % barra de cores
set(get(h,'label'),'string','Evolução do gene'); % label da barra de cores

% A principal diferença é que a concentração parece não atingir o
% equilíbrio tão cedo como o anterior