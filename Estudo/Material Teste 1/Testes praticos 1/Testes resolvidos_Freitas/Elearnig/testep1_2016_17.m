%% 1o Teste Prático de Avaliação Discreta Física Computacional — 2016/2017 24 de março de 2017
%% 1 a)
clear all
close all
clc

L = 67.0; % m
omega = 7.292E-5; % s^-1
phi = 48.865; % graus  em graus usa-se sind(phi)
g = 9.8; % m/s^2
x0 = 2.00; % m
y0 = 0.00; % m
v0 = 0; % m/s
t0 = 0; % s
tf = 200; % s
h = 0.01; % s

t = t0:h:tf;
N = length(t);
x = nan(N,1);
y = nan(N,1);
vx = nan(N,1);
vy = nan(N,1);
x(1) = x0; y(1) = y0;
vx(1) = v0; vy(1) = v0;

% Metódo de Euler-Cromer

for k = 1:(N-1)
    vx(k+1) = vx(k) + (2 * omega * sind(phi) * vy(k) - g / L * x(k)) * h;
    x(k+1) = x(k) + vx(k+1) * h;
    
    vy(k+1) = vy(k) + (-2 * omega * sind(phi) * vx(k) - g / L * y(k)) * h;
    y(k+1) = y(k) + vy(k+1) * h;
end

% % com funcoes anonimas
% fx = @(vx) vx;
% fy = @(vy) vy;
% fvx = @ (x, vy) 2 * omega * sphi * vy - g/L * x;
% fvy = @ (y, vx) -2 * omega * sphi * vx - g/L * y;
% 
% for i = 1:n-1
%     vx(i+1) = vx(i) + fvx(x(i), vy(i)) * h;
%     x(i+1) = x(i) + vx(i+1) * h;
%     
%     vy(i+1) = vy(i) + fvy( y(i), vx(i)) * h;
%     y(i+1) = y(i) + vy(i+1) * h;
% end

plot(x,y,'b.-'), xlabel('x(m)'), ylabel('y(m)'), title('Pêndulo de Foucalt')


%% 1 b)
clear all
close all
clc

L = 67.0; % m
omega = 7.292E-5; % s^-1
phi = 48.865; % graus  em graus usa-se sind(phi)
g = 9.8; % m/s^2
x0 = 2.00; % m
y0 = 0.00; % m
v0 = 0; % m/s
t0 = 0; % s
tf = 500; % s
h = 0.01; % s

t = t0:h:tf;
N = length(t);
x = nan(N,1);
y = nan(N,1);
vx = nan(N,1);
vy = nan(N,1);
ang = nan(N,1); ang(1) = 0;
x(1) = x0; y(1) = y0;
vx(1) = v0; vy(1) = v0;

% Metódo de Euler-Cromer

fx = @(vx) vx;
fy = @(vy) vy;
fvx = @ (x, vy) 2 * omega * sind(phi) * vy - g/L * x;
fvy = @ (y, vx) -2 * omega * sind(phi) * vx - g/L * y;

ind_max = [];       %matriz que vai conter os indices dos maximos locais de x

for i = 1:N-1
    vx(i+1) = vx(i) + fvx(x(i), vy(i)) * h;
    x(i+1) = x(i) + vx(i+1) * h;
    
    vy(i+1) = vy(i) + fvy( y(i), vx(i)) * h;
    y(i+1) = y(i) + vy(i+1) * h;
    if i > 1
        if x(i) > x(i-1) && x(i) > x(i+1)       %se x for maior que o x anterior e o x seguinte
            ind_max = [ind_max i];
        end
    end
end

nmax = length(ind_max);     %numero de maximos locais
t_max = zeros(1, nmax);
theta_max = zeros(1, nmax);

for j = 1:length(ind_max)
    ind = ind_max(j);           %indice atual
    t_max(j) = t(ind);
    theta_max(j) = mod (atan2 (y(ind), x(ind)) , 2*pi);     %calculo de theta usando mod
end

plot(t_max, theta_max, '*k'), xlabel('tempo (s)'), ylabel('theta (rad)'), title('Máximos locais de x')
lsline
p = polyfit(t_max, theta_max, 1);
declive = p(1);
ordenada_na_origem = p(2);

T_precessao = 2*pi / abs(declive);      %periodo de precessão em segundos

disp(['periodo: ',num2str(T_precessao/3600),' horas'])

%% 2 a)
clear all
close all
clc

L = 0.100;
ro = 2.70*10^3;
sigma = 5.6703*10^-8;
A = L^2;
m = A*L / ro;
c = 0.91*10^3;
Tc = 0;

ti = 0;
tf = 100;
h = 0.01;
t = ti:h:tf;
N = length(t);

T = zeros(1, N);

T(1) = 2000;

% metodo de crank-nicolson
for i = 1:N-1
    %vetor com os coeficientes do polinomio de quarto grau
    C = [ sigma * A/ (m * c) * h/2 , 0, 0, 1, - (T(i) - sigma * A/ (m * c) * h/2 *(T(i)^4 - 2 * Tc^4))];
    sol = roots(C);
    [valor, indice] = min(abs(sol-T(i)));
    T(i+1) = sol(indice);
end

plot( t, T, 'b')
xlabel('tempo(s)')
ylabel('Temperatura(K)')
title('Arrefecimento de um bloco de alumínio')

%% 2) b)
clear all
close all
clc

L = 0.100;
ro = 2.70*10^3;
sigma = 5.6703*10^-8;
A = L^2;
m = A*L / ro;
c = 0.91*10^3;
Tc = 0;
T0 = 2000;

ti = 0;
tf = 1000;
hs = [ 0.01 0.02 0.025 0.04 0.05];
n = length(hs);

Tf_fin = (T0^-3 + 3*sigma*A / (m*c) * tf)^(-1/3);            %temperatura final analitica

erro = zeros(1, n);            %erro para cada h       

for j = 1:n
    h = hs(j);
    t = ti:h:tf;
    N = length(t);
    T = zeros(1, N);
    T(1) = T0;
    
    for i = 1:N-1
        %vetor com os coeficientes do polinomio de quarto grau
        C = [ sigma * A/ (m * c) * h/2 , 0, 0, 1, - (T(i) - sigma * A/ (m * c) * h/2 *(T(i)^4 - Tc^4))];
        sol = roots(C);
        [valor, indice] = min(abs(sol-T(i)));
        T(i+1) = sol(indice);
    end
    
    erro(j) = abs(Tf_fin - T(end));  %modulo da diferença entre o Tf obtido e o Tf analitico
    
end

plot( log(hs), log(erro), '*b')
lsline
xlabel('log(h) ')
ylabel('log(Erro)')
title('Arrefecimento de um bloco de alumínio')

p = polyfit(log(hs), log(erro), 1);
ordem_do_metodo = p(1);         %ordem do metodo é o declive da reta
disp(['Ordem do metódo: ',num2str(round(ordem_do_metodo)),' '])

%% 2 c)
clear all
close all
clc

L = 0.1; % m
A = L*L;
p = 2.70E3; % kg/m^3
m = p * A * L; % kg
sigma = 5.6706E-8; % W m^-2 K^-4
c = 0.91E3; %J kg^-1 K^-1
T0 = 310; % K


t0 = 0; % s
tf = 2*3600; % 2 horas em s
h = 0.01; %s

t = t0:h:tf; %s
N = length(t);
T = nan(N,1);
T(1) = T0;
% Tc = nan(N,1);
% Tc(1) = 283 + 1.0*10^-3 * 0;


fT = @(t, T) -sigma*A/(m*c) * (T^4 - (283 + 1.0*10^-3 * t)^4);
    
for i = 1:N-1
    
    r1T = fT(t(i), T(i));
    
    r2T = fT(t(i) + h/2, T(i) + r1T * h/2);
    
    r3T = fT(t(i) + h/2, T(i) + r2T * h/2);
    
    r4T = fT(t(i) + h, T(i) + r3T * h);

    
    T(i+1) = T(i) + 1/6*(r1T + 2*r2T + 2*r3T + r4T)*h;

end

figure(2)
plot(t,T,'b.-'),xlabel('t(s)'), ylabel('T(k)'), title('Temperatura em função do tempo nas primeiras duas horas')

%% 2 c) outra forma
clear all
close all
clc

L = 0.1; % m
A = L*L;
p = 2.70E3; % kg/m^3
m = p * A * L; % kg
sigma = 5.6706E-8; % W m^-2 K^-4
c = 0.91E3; %J kg^-1 K^-1
T0 = 310; % K


t0 = 0; % s
tf = 2*3600; % 2 horas em s
h = 0.01; %s

t = t0:h:tf; %s
N = length(t);
T = nan(N,1);
T(1) = T0;


fT = @(T, Tc) -sigma*A/(m*c) * (T^4 - Tc^4);
Tc = @(t) 283 + 1.0*10^-3 * t;
    
for i = 1:N-1
    
    r1T = fT(T(i), Tc(t(i)));
    r2T = fT(T(i) + r1T * h/2, Tc(t(i)+0.5*h));
    r3T = fT(T(i) + r2T * h/2, Tc(t(i))+0.5*h);
    r4T = fT(T(i) + r3T * h, Tc(t(i))+h);
    T(i+1) = T(i) + 1/6*(r1T + 2*r2T + 2*r3T + r4T)*h;
    
end

plot(t,T,'b.-'),xlabel('t(s)'), ylabel('T(k)'), title('Temperatura em função do tempo nas primeiras duas horas')
