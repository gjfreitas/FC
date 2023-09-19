%% 2o Teste Prático de Avaliação Discreta Física Computacional — 2017/2018
%% 2 a) e b)
clear all
close all
clc

% Constantes
L = 2.00;
c = 1.50;
a = 1.00 * 10^(-3);
C = 1;

x0 = 0;
xf = L;
t0 = 0;
tf = 15;

dx = 0.1;
dt = C*dx/c;

x = x0:dx:xf;
t = t0:dt:tf;

Nx = length(x);
Nt = length(t);

y = zeros(Nx,Nt);
C = c*dt/dx;

% Condições iniciais
y(:,1) = 0;
pd = a*x.*(L-x);

% Condições fronteira
y(1,:) = 0;
y(Nx,:) = 0;

f = @(x,t) a*(2*c^2-x*(L-x))*sin(t);

for n = 1: Nt-1
    for i = 2:Nx-1
        if n == 1
            y(i,2) = y(i,1) + (C^2/2)*(y(i-1,1)-2*y(i,1)+y(i+1,1)) + pd(i)*dt +(1/2)*dt^2*f(x(i),t(1));
        else
            y(i,n+1) = - y(i,n-1) + 2*y(i,n) + C^2*(y(i-1,n)-2*y(i,n)+y(i+1,n)) + dt^2*f(x(i),t(n));
        end
    end
end

figure(1)
mesh(t,x,y), xlabel('t'), ylabel('x'), zlabel('y')
figure(2)
contourf(x,t,y'), xlabel('x'), ylabel('t'), zlabel('y')
colorbar; % barra de cores


%% 2 c)
clear all
close all
clc

% Constantes
L = 2.00;
c = 1.50;
a = 1.00 * 10^(-3);
C = 1;

x0 = 0;
xf = L;
t0 = 0;
tf = 7.2;

dx = 0.2;
dt = C*dx/c;

x = x0:dx:xf;
t = t0:dt:tf;

Nx = length(x);
Nt = length(t);

y = zeros(Nx,Nt);
C = c*dt/dx;

% Condições iniciais
y(:,1) = 0;
pd = a*x.*(L-x);

% Condições fronteira
y(1,:) = 0;
y(Nx,:) = 0;

f = @(x,t) a*(2*c^2-x*(L-x))*sin(t);

for n = 1: Nt-1
    for i = 2:Nx-1
        if n == 1
            y(i,2) = y(i,1) + (C^2/2)*(y(i-1,1)-2*y(i,1)+y(i+1,1)) + pd(i)*dt +(1/2)*dt^2*f(x(i),t(1));
        else
            y(i,n+1) = - y(i,n-1) + 2*y(i,n) + C^2*(y(i-1,n)-2*y(i,n)+y(i+1,n)) + dt^2*f(x(i),t(n));
        end
    end
end

figure(1)
mesh(t,x,y), xlabel('t'), ylabel('x'), zlabel('y')
figure(2)
contourf(x,t,y'), xlabel('x'), ylabel('t'), zlabel('y')
colorbar; % barra de cores

y_analitico = a*x.*(L-x)*sin(7.2);
figure(3)
plot(x,y(:,Nt),'o',x,y_analitico,'ro-');
legend('Y obtido','Y analitico')

erro = max(abs(y(:,Nt)'- y_analitico));
disp(['Erro: ',num2str(erro)])

%% 2 d)
clear all
close all
clc

% Constantes do sistema
L = 2.00;
c = 1.50;
a = 1.00 * 10^(-3);
C = 1;

% Pré-alocações
x0 = 0;
xf = L;
t0 = 0;
tf = 7.2;

dx_vec = [0.400 0.200 0.100 0.050 0.025];

NN = length(dx_vec);
erro_vec = zeros(1,NN);

dx_max = max(dx_vec);
dx_min = min(dx_vec);

Nx_vec = (xf-x0)./dx_vec + 1;

Nmin = (xf-x0)/dx_max + 1;
Nmax = (xf-x0)/dx_min + 1;

dt_vec = C*dx_vec./c;
dt_min = min(dt_vec);
Nt_vec = (tf-t0)./dt_vec +1;

x = x0:dx_min:xf; % guardar toda a memória utilizada antes do for
t = t0:dt_min:tf; % guardar toda a memória utilizada antes do for


Nx_max = length(x); % guardar toda a memória utilizada antes do for
Nt_max = length(t); % guardar toda a memória utilizada antes do for

y = zeros(Nx_max,Nt_max); % guardar toda a memória utilizada antes do for

f = @(xi,ti) a*(2*c^2-xi*(L-xi))*sin(ti);


for idx = 1:NN
    dx = dx_vec(idx);
    dt = dt_vec(idx);
    Nx = Nx_vec(idx);
    Nt = Nt_vec(idx);
    
    x = x0:dx:xf;
    t = t0:dt:tf;
    
    % Solução analítica
    y_analitico = a*x.*(L-x)*sin(7.2);
    
    % Condições iniciais
    y(1:Nx,1) = 0;
    pd = a*x.*(L-x);
    
    % Condições fronteira
    y(1,1:Nt) = 0;
    y(Nx,1:Nt) = 0;

    % Resolução numérica
    for n = 1: Nt-1
        for i = 2:Nx-1
            if n ==1
                y(i,2) = y(i,1) + (C^2/2)*(y(i-1,1)-2*y(i,1)+y(i+1,1)) + pd(i)*dt +(1/2)*dt^2*f(x(i),t(1));
            else
                y(i,n+1) = - y(i,n-1) + 2*y(i,n) + C^2*(y(i-1,n)-2*y(i,n)+y(i+1,n)) + dt^2*f(x(i),t(n));
            end
        end
    end
    
    erro_vec(idx) = max(abs(y_analitico(1:Nx)-y(1:Nx,Nt)'));
end

plot(log(dx_vec),log(erro_vec));
p = polyfit(log(dx_vec),log(erro_vec),1);

disp(['O log do erro é proporcional ao log do dx, sendo o declive da reta m = ',num2str(p(1)),' . Logo, deve existir uma relação: erro = const*dx^(',num2str(p(1)),').']) 