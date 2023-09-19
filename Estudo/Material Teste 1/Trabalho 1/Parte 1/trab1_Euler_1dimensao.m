%% Trabalho 1 - Movimento a uma dimensão
%% A - Queda de uma pedra 
clc
clear all
close all

% f= ma <=> -mg = ma <=> (d^2/dt^2)(z) = -g     a = (d^2/dt^2)(z) uma
% dimensão(z)

% a) (d^2/dt^2)(z) = -g 

m = 150*10^-3; % kg
g = 9.8; % m/s^2
alt = 3; % m
v0 = 0; % m/s a pedra é largada
z0 = 2*alt; % m altura incial (segundo andar)

% b) Não, pois a aceleração é um equação de segunda ordem e o Metódo
% resolve equações de primeira ordem. Logo, não a podemos usar diretamente.

% c) Transforma-mos um equação linear de segunda ordem num sistema de duas
% equações lineares de primeira ordem
% dv/dt = -g e dz/dt = v  o que leva a:
% v(k+1) = v(k) +(-g)*h <- Metódo de Euler
% z(k+1) = z(k) + v(k)*h <- Metódo de Euler

t0 = 0; %s
tf = 1.5; %s
h = 0.2; %s

t = t0:h:tf; %s
N = length(t);
z = nan(1,N);
v = nan(1,N);
z(1) = z0; v(1) = v0;

for k = 1:N-1
    v(k+1) = v(k) - g*h;
    z(k+1) = z(k) + v(k) * h;
end

% e)

vz = -g.*t+v0; % solução analitica

figure(2)
plot(t,v,'.r-',t,vz,'ko--')
title("velocidade em função do tempo")
xlabel("tempo(s)")
ylabel("velocidade(m/s^2)")

% g)
zz = z0 - 0.5 * g .* t.^2;
figure(2)
plot(t,z,'.r-',t,zz,'ko--')
title("posição em função do tempo")
xlabel("tempo(s)")
ylabel("altura(m)")

%% f) (Estimativa numérica do instante em que a pedra cai no chão e da sua velocidade.)
clc
clear all
close all

m = 150*10^-3; % kg
g = 9.8; % m/s^2
alt = 3; % m
v0 = 0; % m/s a pedra é largada
z0 = 2*alt; % m altura incial (segundo andar)
t0 = 0; %s
tf = 1.5; %s
h = 0.2; %s

t = t0:h:tf; %s
N = length(t);
z = nan(1,N);
v = nan(1,N);
z(1) = z0; v(1) = v0;

for k = 1:N-1
    v(k+1) = v(k) - g*h;
    z(k+1) = z(k) + v(k) * h;
    
    if z(k) < 0
        break % penultimo valor de z positivo e ultimo negativo, isto faz-se para se poder usar interp1
    end
end

t_impacto = interp1(z(end:-1:end-1),t(end:-1:end-1),0)
v_final = v(end)


%% B - Oscilador harmónico simples (aceleração variavel)

% a) F = -kx e F = ma <=> ma = -kx <=> a = (-k/m)x <=> (d^2/dt^2)(x) = (-k/m)x

% b)Para usar o metódo de Euler transforma-mos um equação linear de segunda ordem num sistema de duas
% equações lineares de primeira ordem
% Temos então : dv/dt = (-k/m)x e dx/dt = v
% Logo as equções para o Metódo de Euler : v(k+1) = v(k) + ((-k/m)x)*h e x(k+1) = x(k) + v(k)*h

clc
clear all 
clc 

K = 2.5; %N/m
m = 0.5; %kg
x0 = 10*10^-2; %m
v0 = 0;

% c)

h = 0.1; %s
% h = 0.01; %s
% Quanto menor o h maior vai ser a precisão dos cálculos (mais proximos dos valores exatos)

t0 = 0; %s
tf = 10; %s
t = t0:h:tf;

N = length(t);
v = zeros(1,N);
v(1) = v0;
x = zeros(1,N);
x(1) = x0;

for k = 1:N-1
    v(k+1) = v(k) + (-K/m)*x(k)*h;
    x(k+1) = x(k) + v(k)*h;
end

% d)
w = sqrt(K/m);
x_exact = x0*cos(w*t); % solução analitica
v_exact = -x0*w*sin(w*t); % solução analitica

figure(2)
plot(t,x,'-b',t,x_exact,'ok--')
title("posição em função do tempo")
xlabel("tempo(s)")
ylabel("posição(m)")

figure(3)
plot(t,v,'-b',t,v_exact,'ok--')
title("velocidade em função do tempo")
xlabel("tempo(s)")
ylabel("velocidade(m/s^2)")


%% e)
K = 2.5; %N/m
m = 0.5; %kg
x0 = 10*10^-2; %m

h = [1E-1, 5E-2,1E-2, 5E-3,1E-3, 5E-4,1E-4]; %s
% h = [1E-2, 5E-3,1E-4, 5E-4,1E-4, 5E-5,1E-5]; %s
erro = zeros(1,length(h));

for i = 1:length(h)
    t0 = 0; %s
    tf = 10; %s
    t = t0:h(i):tf;

    N = length(t);
    v = zeros(1,N);
    v(1) = v0;
    x = zeros(1,N);
    x(1) = x0;
    
    for k = 1:N-1
        v(k+1) = v(k) + (-K/m)*x(k)*h(i);
        x(k+1) = x(k) + v(k)*h(i);
    end
    
    erro(i) = abs(x0-max(x));
    
end

plot(log10(h),log10(erro),'ko')
lsline; %traçar uma recta sobre a curva que obteve
aux = polyfit(log10(h),log10(erro),1);  %1 corresponde a ordem do polinomio, como é uma reta a ordem = 1
declive = aux(1);
disp(['Declive : ',num2str(declive),' '])

% Se usarmos valores de h mais pequenos a reta fica mais perto dos valores de erro e o valor do declive aproxima-se de 1
    


