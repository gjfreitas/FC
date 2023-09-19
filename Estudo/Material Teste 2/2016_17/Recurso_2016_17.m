%% Exame Prático de Recurso Física Computacional — 2016/2017
%% 2 a)
clear all
close all
clc

l = 0.0027;
alfa = pi/4;
a = 0.0025;
tmax = pi/2-alfa;
t0 = 0;
h = 0.0054;
t = 0:h:tmax; % t = theta

S0 = 0;
X0 = 0;
SF = 1;
dS0 = 2*l^2/(a*h);
dX0 = 0;

N = length(t);
X = nan(1,N); X(1) = 0;
S = X;
dX = X;
dS = nan(1,N); dS(1) = dS0;

Q = nan(1,N); W = nan(1,N);
Q(1) = a/(l^2) * (X0 + h/a) - sin(t0)/(a); % aqui tirei o x S0 pq dividir por 0 é impossivel
W(1) = (a/l)^2 * (X0 + h/a) - sin(t0)/Q(1);


for k = 1:N-1
    
    dS(k+1) = dS(k) + (cos(t(k))/W(k)) * h;
    S(k+1) = S(k) + dS(k+1)*h;
    
    dX(k+1) = dX(k) + (sin(t(k))/W(k)) * h;
    X(k+1) = X(k) + dX(k+1) * h;
    
    Q(k+1) = a/(l^2) * (X(k+1) + h/a) - sin(t(k+1))/(a*S(k+1));
    W(k+1) = (a/l)^2 * (X(k+1) + h/a) - sin(t(k+1))/Q(k+1);
    
end

disp(['S(theta_max) = ',num2str(S(end))])

%% 2 b)

clear all
close all
clc

l = 0.0027;
alfa = pi/4;
a = 0.0025;
tmax = pi/2-alfa;
t0 = 0;

S0 = 0;
X0 = 0;
dX0 = 0;


% shooting
B = 1;
guess = [0.0020 0.0045];
tol = 1E-12;
nshots = 1000;

for j = 1:nshots

    h = guess(j);
    t = 0:h:tmax; % t = theta
    
    dS0 = 2*l^2/(a*h);

    N = length(t);
    X = nan(1,N); X(1) = 0;
    S = X;
    dX = X;
    dS = nan(1,N); dS(1) = dS0;
    Q = nan(1,N); W = nan(1,N);
    Q(1) = a/(l^2) * (X0 + h/a) - sin(t0)/(a); % aqui tirei o x S0 pq dividir por 0 é impossivel
    W(1) = (a/l)^2 * (X0 + h/a) - sin(t0)/Q(1);

    for k = 1:N-1
        dS(k+1) = dS(k) + (cos(t(k))/W(k)) * h;
        S(k+1) = S(k) + dS(k+1)*h;

        dX(k+1) = dX(k) + (sin(t(k))/W(k)) * h;
        X(k+1) = X(k) + dX(k+1) * h;

        Q(k+1) = a/(l^2) * (X(k+1) + h/a) - sin(t(k+1))/(a*S(k+1));
        W(k+1) = (a/l)^2 * (X(k+1) + h/a) - sin(t(k+1))/Q(k+1);
    end
    
    result(j) = S(end);
    diff = B-result(j);
       
    if j>= 2
        m = (result(j)-result(j-1))/(guess(j)-guess(j-1));
        guess(j+1) = guess(j)+(diff)/m;
        if abs(guess(j) - guess(j-1)) < tol
            break
        end
    end
end

sol = result(j);
hc = guess(j);
disp(['S(theta_max) = ',num2str(sol)])
disp(['h = ',num2str(hc)])

%% 3 a)
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
colorbar

%% 3 b)
clear all
close all
clc

L = 2;
D = 1;
C = 1;
L2 = 4;
dx = 0.1;
x = -L/2:dx:L/2;
x2 = -L2/2:dx:L2/2;
dt = 0.001;
t = 0:dt:1;

Nx = length(x);
Nx2 = length(x2);
Nt = length(t);


p = nan(Nx,Nt);
p2 = nan(Nx2,Nt);

% Condiçoes iniciais
p(:,1) = 2*cos(pi*x/L);
p2(:,1) = 2*cos(pi*x2/L2);

% Condicoes de fronteira
p(1,:) = 0;
p(end,:) = 0;
p2(1,:) = 0;
p2(end,:) = 0;

for j = 1:Nt-1
    for i = 2:Nx-1
        p(i,j+1) = p(i,j) + dt*(D/(dx^2) * (- 2*p(i,j)+ p(i-1,j) + p(i+1,j))+C*p(i,j)) ;
        p2(i,j+1) = p2(i,j) + dt*(D/(dx^2) * (- 2*p2(i,j)+ p2(i-1,j) + p2(i+1,j))+C*p2(i,j)) ;
    end
end

% L = 2
figure(1)
mesh(t,x,p)
xlabel('t')
ylabel('x')
zlabel('p')
title('L = 2')
figure(2)
contourf(x,t,p')
colorbar
title('L = 2')

% L = 4
figure(3)
mesh(t,x2,p2)
xlabel('t')
ylabel('x')
zlabel('p')
title('L = 4')
figure(4)
contourf(x2,t,p2')
colorbar
title('L = 4')

disp('Pelo grafico, podemos ver que quando L = 4, o comportamento se alterou')