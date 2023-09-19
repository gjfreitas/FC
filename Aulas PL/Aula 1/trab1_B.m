%% Trabalho 1 - Movimento a uma dimensão
%% B - Oscilador harmónico simples (aceleração variavel)

% a) F = -kx e F = ma <=> ma = -kx <=> a = (-k/m)x <=> (d^2/dt^2)(x) = (-k/m)x

% b)Para usar o metódo de Euler transforma-mos um equação linear de segunda ordem num sistema de duas
% equações lineares de primeira ordem
% Temos então : dv/dt = (-k/m)x e dx/dt = v (Importante!!!!!)
% Logo as equções para o Metódo de Euler : v(k+1) = v(k) + ((-k/m)x)*h e x(k+1) = x(k) + v(k)*h  (Importante!!!!!)

%% c) Parte 1 - h = 0.1
clc
clear all 
clc 
K = 2.5; %N/m
m = 0.5; %kg
x0 = 10*10^-2; %m
t0 = 0; %s
h = 0.1; %s
tf = 10; %S
t = t0:h:tf;
v0 = 0;

N = length(t);
v_euler = zeros(1,N);
v_euler(1) = v0;
x_euler = zeros(1,N);
x_euler(1) = x0;

for k = 1:(N-1) %N-1 para os comprimentos de v e t serem iguais
    v_euler(k+1) = v_euler(k)+((-K/m)*x_euler(k))*h; %((-K/m)*x)*h pois a derivada da velocidade é (-K/m)*x
    x_euler(k+1) = x_euler(k)+v_euler(k)*h; %v(k)*h pois a derivada de x é v
end

% figure(1)
% plot(t,x_euler,'-b')
% title("posição em função do tempo")
% xlabel("tempo(s)")
% ylabel("posição(m)")

% d) parte analitica
% w^2 = K/m
w = sqrt(K/m);
x_exact = x0*cos(w.*t);
v_exact = -x0.*w.*sin(w.*t);

figure(2)
plot(t,x_euler,'-b',t,x_exact,'ok--')
title("posição em função do tempo")
xlabel("tempo(s)")
ylabel("posição(m)")

figure(3)
plot(t,v_euler,'-b',t,v_exact,'ok--')
title("velocidade em função do tempo")
xlabel("tempo(s)")
ylabel("velocidade(m/s^2)")

%% c) Parte 1 - h = 0.01
clc
clear all 
clc 
K = 2.5; %N/m
m = 0.5; %kg
x0 = 10*10^-2; %m
t0 = 0; %s
h = 0.01; %s
tf = 10; %S
t = t0:h:tf;
v0 = 0;

N = length(t);
v_euler = zeros(1,N);
v_euler(1) = v0;
x_euler = zeros(1,N);
x_euler(1) = x0;

for k = 1:(N-1) %N-1 para os comprimentos de v e t serem iguais
    v_euler(k+1) = v_euler(k)+((-K/m)*x_euler(k))*h; %((-K/m)*x)*h pois a derivada da velocidade é (-K/m)*x
    x_euler(k+1) = x_euler(k)+v_euler(k)*h; %v(k)*h pois a derivada de x é v
end

% w^2 = K/m
w = sqrt(K/m);

x_exact = x0*cos(w.*t);
v_exact = -x0.*w.*sin(w.*t);

figure(1)
subplot(1,2,1)
plot(t,x_euler,'-b',t,x_exact,'k--')
title("posição em função do tempo")
xlabel("tempo(s)")
ylabel("posição(m)")
subplot(1,2,2)
plot(t,v_euler,'-b',t,v_exact,'k--')
title("velocidade em função do tempo")
xlabel("tempo(s)")
ylabel("velocidade(m/s^2)")

% Logo, quanto menor for o incremento (h) entre os tempos menor vai ser a
% diferença entre os valores exactos e os valores utilizando o Metódo de Euler

%% e) 
clc
clear all 
clc 
K = 2.5; %N/m
m = 0.5; %kg
x0 = 10*10^-2; %m
t0 = 0; %s
tf = 10; %S

h = [1E-1, 5E-2,1E-2, 5E-3,1E-3, 5E-4,1E-4]; %s
% h = [1E-2, 5E-3,1E-4, 5E-4,1E-4, 5E-5,1E-5]; %s
erro = zeros(1,length(h));

for index_h = 1:length(h)
    t = t0:h(index_h):tf;
    v0 = 0;

    N = length(t);
    v_euler = zeros(1,N);
    v_euler(1) = v0;
    x_euler = zeros(1,N);
    x_euler(1) = x0;

    for k = 1:(N-1) %N-1 para os comprimentos de v e t serem iguais
        v_euler(k+1) = v_euler(k)+((-K/m)*x_euler(k))*h(index_h); %((-K/m)*x)*h pois a derivada da velocidade é (-K/m)*x
        x_euler(k+1) = x_euler(k)+v_euler(k)*h(index_h); %v(k)*h pois a derivada de x é v
    end
    
    erro(index_h) = abs(x0-max(x_euler));
    
end

plot(log10(h),log10(erro),'ko')
lsline; %traçar uma recta sobre a curva que obteve
aux = polyfit(log10(h),log10(erro),1)  %1 corresponde a ordem do polinomio, como é uma reta a ordem = 1
declive = aux(1) 

% Se usarmos valores de h mais pequenos a reta fica mais perto dos valores de erro e o valor do declive aproxima-se de 1


