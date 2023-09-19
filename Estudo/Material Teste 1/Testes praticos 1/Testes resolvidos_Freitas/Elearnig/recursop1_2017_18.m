%% Exame Prático de Recurso — Parte 1 Física Computacional — 2017/2018 3 de julho de 2018

% 1 a)
clear all
close all
clc

% Constantes
m = 0.5; %kg
I = 10^-4; %kg m^2
K = 5; %N/m
lambda = 10^-3; %N m
eps = 10^-2; %N

%Condições iniciais
z0 = 0.1; %m
v0 = 0; %m/s
theta0 = 0; %rads
w0 = 0; %rads/s

h = 0.01; %s
t = 0:h:100; %s

N = length(t);
z = nan(1,N); z(1) = z0;
v = nan(1,N); v(1) = v0;
w = nan(1,N); w(1) = w0;
theta = nan(1,N); theta(1) = theta0;

A_CN = [ h*K/(2*m), 1, h*eps/(4*m), 0; 1, -h/2, 0, 0; h*eps/(4*I), 0, h*lambda/(2*I), 1; 0, 0, 1, -h/2];
b_CN = [v0 - h*K/(2*m) * z0 - h*eps/(4*m) * theta0; z0 + h/2 * v0; w0 - h*eps/(4*I) * z0 - h*lambda/(2*I) * theta0; theta0 + h/2 * w0];

for k = 1:(N-1)
    J_CN = linsolve(A_CN,b_CN);
    z(k+1) = J_CN(1);
    v(k+1) = J_CN(2);
    theta(k+1) = J_CN(3);
    w(k+1) = J_CN(4);
    
    b_CN = [J_CN(2) - h*K/(2*m) * J_CN(1) - h*eps/(4*m) * J_CN(3); J_CN(1) + h/2 * J_CN(2); J_CN(4) - h*eps/(4*I) * J_CN(1) - h*lambda/(2*I) * J_CN(3); J_CN(3) + h/2 * J_CN(4)];
end

figure(1)
subplot(2,2,1)
plot(t,z), xlabel('t(s)'), ylabel('z')
subplot(2,2,2)
plot(t,v), xlabel('t(s)'), ylabel('v')
subplot(2,2,3)
plot(t,theta), xlabel('t(s)'), ylabel('theta')
subplot(2,2,4)
plot(t,w), xlabel('t(s)'), ylabel('w')



count = 0;

for k = 2:(N-1) % começamos em 2 pois se começarmos em 1 não temos nenhum anterior para comparar (em matlab não ha posição 0)
    if z(k-1) <= z(k) && z(k) >= z(k+1)
        count = count + 1;
        
        aux = lagr(t(k-1:k+1),z(k-1:k+1));
        t_max(count) = aux(1);

    end
end

p = mean(diff(t_max)) ; %periodo
disp(['Período médio de Frequência: ',num2str(p),' s'])
freq = 1/p;
disp(['Frequência =  ',num2str(freq),' Hz '])


% c)
figure(2)
E = 0.5* K .* z.^2 + 0.5 * m .* v .^2;
plot(t,E), xlabel('t(s)'), ylabel('E(J)')


count = 0;
for k = 2:(N-1) % começamos em 2 pois se começarmos em 1 não temos nenhum anterior para comparar (em matlab não ha posição 0)
    if E(k-1) <= E(k) && E(k) >= E(k+1)
        count = count + 1;
        
        aux = lagr(t(k-1:k+1),E(k-1:k+1));
        t_max(count) = aux(1);
    end
end
p = mean(diff(t_max)) ; %periodo
disp(['Período médio de Energia: ',num2str(p),' s'])
disp(['Intervalos de tempo fixos em q o sistema muda de oscilações puramente longitudinais para oscilações puramente torsionais : ',num2str(p/2),' s'])


%% 2 a)
clear all
close all
clc

a = 2;
b = 0.74;
c = 0.5;
x0 = 0.8;
y0 = 0.3;
h = 0.01;
tf = 100;
t = 0:h:tf;
N = length(t);

x = nan(1,N); x(1) = x0;
y = nan(1,N); y(1) = y0;
vx = nan(1,N); vx(1) = 0;
vy = nan(1,N); vy(1) = 0;

fx = @(x,y) x*(1-x) - (a*x*y)/(x+y);
fy = @(x,y) (b*x*y)/(x+y) - c*y;

for k = 1:N-1

    r1x = fx(x(k),y(k));
    r1y = fy(x(k),y(k));
    
    r2x = fx(x(k)+r1x*h/2, y(k)+r1y*h/2);
    r2y = fy(x(k)+r1x*h/2, y(k)+r1y*h/2);
        
    r3x = fx(x(k)+r2x*h/2, y(k)+r2y*h/2);
    r3y = fy(x(k)+r2x*h/2, y(k)+r2y*h/2);
    

    r4x = fx(x(k)+r3x*h/2, y(k)+r3y*h);
    r4y = fy(x(k)+r3x*h/2, y(k)+r3y*h);
       
    x(k+1) = x(k) + 1/6 *(r1x + 2*r2x + 2*r3x +r4x)*h;
    y(k+1) = y(k) + 1/6 *(r1y + 2*r2y + 2*r3y +r4y)*h;

    
end
figure(1)
plot(t,x), xlabel('t'), ylabel('x')
figure(2)
plot(t,y), xlabel('t'), ylabel('y')

%% 2 b)

x0 = 0.8;
y0 = 0.3;
h = 0.01;
tf = 100;
t0 = 0;
t = 0:h:tf;


% ode45
reltol = 3E-14;
abstol_1 = 1E-13;
abstol_2 = 1E-13;
options = odeset('RelTol',reltol,'AbsTol',[abstol_1 abstol_2]);
[t_ode45,sol] = ode45(@f,[t0 tf],[x0 y0],options);

figure(1)
plot(t_ode45,sol(:,1),'r--'), xlabel('t'), ylabel('x')
figure(2)
plot(t_ode45,sol(:,2),'r--'), xlabel('t'), ylabel('y')

function derivadas = f(t,solucao)
a = 2;
b = 0.74;
c = 0.5;
derivadas = zeros(2,1);

% O vetor solucao tem os valores de x e v para o tempo t em que a função é chamada pela rotina ode45.
derivadas(1) = solucao(1) * ( 1 - solucao(1)) - a * solucao(1) * solucao(2) / (solucao(1) + solucao(2));
% Escreva acima a expressão da derivada de x em função de solucao(1) e de solucao(2).

derivadas(2) = b * solucao(1) * solucao(2) / (solucao(1) + solucao(2)) - c * solucao(2);
% Escreva a acima expressão da derivada de v em função de solucao(1) e de solucao(2). 
% Se achar que torna o programa mais claro, pode incluir na função as linhas,
% x = solucao(1);
% y = solucao(2);
end