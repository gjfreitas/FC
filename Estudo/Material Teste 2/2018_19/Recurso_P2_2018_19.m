%% Exame Prático de Recurso — Parte 2 Física Computacional — 2018/2019
%% 1 a)
clear all
close all
clc

L = 2;
D = 1;
dx = 0.1;
x = -L/2:dx:L/2;
dt = 0.001;
t = 0:dt:1;

Nx = length(x);
Nt = length(t);

p = nan(Nx,Nt);

% Condiçoes iniciais
p(:,1) = 2*cos(pi*x/L);

% Condicoes de fronteira
p(1,:) = 0;
p(end,:) = 0;

for j = 1:Nt-1
    for i = 2:Nx-1
        p(i,j+1) = p(i,j) + D*dt/(dx^2) * (- 2*p(i,j)+ p(i-1,j) + p(i+1,j)); % slide 4 - Aula 6 (2020/21)
    end
end

figure(1)
mesh(t,x,p), xlabel('t'), ylabel('x'), zlabel('p')
figure(2)
contourf(x,t,p'), xlabel('x'), ylabel('t'), zlabel('p')

%% 1 b)
clear all
close all
clc

L = 2;
D = 1;
C = 1;
C2 = 3;
dx = 0.1;
x = -L/2:dx:L/2;
dt = 0.001;
t = 0:dt:1;

Nx = length(x);
Nt = length(t);

p = nan(Nx,Nt);

% Condiçoes iniciais
p(:,1) = 2*cos(pi*x/L);

% Condicoes de fronteira
p(1,:) = 0;
p(end,:) = 0;
p2 = p;

for j = 1:Nt-1
    for i = 2:Nx-1
        p(i,j+1) = p(i,j) + dt*(D/(dx^2) * (- 2*p(i,j)+ p(i-1,j) + p(i+1,j))+C*p(i,j)) ;
        p2(i,j+1) = p2(i,j) + dt*(D/(dx^2) * (- 2*p2(i,j)+ p2(i-1,j) + p2(i+1,j))+C2*p2(i,j)) ;
    end
end

% C = 1
figure(1)
mesh(t,x,p), xlabel('t'), ylabel('x'), zlabel('p'), title('C = 1')
figure(2)
contourf(x,t,p'), xlabel('x'), ylabel('t'), zlabel('p'), title('C = 1')


% C = 3
figure(3)
mesh(t,x,p2), xlabel('t'), ylabel('x'), zlabel('p'), title('C = 3')
figure(4)
contourf(x,t,p2'), xlabel('x'), ylabel('t'), zlabel('p'), title('C = 3')
colorbar

disp('Pelo grafico, podemos ver que quando C = 3, o comportamento se alterou')

%% 1 c)
clear all
close all
clc

L = 2;
D = 1;
C = 1;
C2 = 2.5;
dx = 0.1;
x = -L/2:dx:L/2;
dt = 0.001;
t = 0:dt:1;

Nx = length(x);
Nt = length(t);

p = nan(Nx,Nt);

% Condiçoes iniciais
p(:,1) = 2*cos(pi*x/L);

% Condicoes de fronteira
p(1,:) = 0;
p(end,:) = 0;
p2 = p;

for j = 1:Nt-1
    for i = 2:Nx-1
        if j == Nt-1
            p(i,j+1) = p(i,j) + dt*(D/(dx^2) * (0)+C*p(i,j)) ;
            p2(i,j+1) = p2(i,j) + dt*(D/(dx^2) * (0)+C2*p2(i,j)) ;
        else
            p(i,j+1) = p(i,j) + dt*(D/(dx^2) * (- 2*p(i,j)+ p(i-1,j) + p(i+1,j))+C*p(i,j)) ;
            p2(i,j+1) = p2(i,j) + dt*(D/(dx^2) * (- 2*p2(i,j)+ p2(i-1,j) + p2(i+1,j))+C2*p2(i,j)) ;
        end
    end
end

% C = 1
figure(1)
mesh(t,x,p), xlabel('x'), ylabel('t'), zlabel('p'), title('C = 1')
figure(2)
contourf(x,t,p'), xlabel('x'), ylabel('t'), zlabel('p'), title('C = 1')
colorbar


% C = 2.5
figure(3)
mesh(t,x,p2), xlabel('x'), ylabel('t'), zlabel('p'), title('C = 2.5')
figure(4)
contourf(x,t,p2'), xlabel('x'), ylabel('t'), zlabel('p'), title('C = 2.5')
colorbar

disp('Pelo grafico, podemos ver que quando C ~ 2.5, o comportamento já se alterou')

%% 2 a) Não temos o ficheiro