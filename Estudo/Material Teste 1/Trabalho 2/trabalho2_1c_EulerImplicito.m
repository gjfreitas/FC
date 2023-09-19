%% Trabalho Prático 2
%% Problema 2.1: Oscilador harmónico simples
clear all
close all 
clc

% c) Metódo de Euler Implicito

K = 1; %N/m
m = 1; %kg
x0 =1; %m
v0 = 0; %m/s

t0 = 0; %s
tf = 50; %s
h = 0.01; %s

t = t0:h:tf; %s
N = length(t);

v = zeros(1,N);
v(1) = v0;

x = zeros(1,N);
x(1) = x0;

% Metódo de Euler Implícito

w = sqrt(K/m);
A = [1,-h; w^2*h,1];
b = [x0;v0];

for k = 1:(N-1)
    z = linsolve(A,b);
    x(k+1) = z(1);
    v(k+1) = z(2);
    
    b = [z(1);z(2)];
end

x_exact = x0*cos(w.*t);
v_exact = -x0.*w.*sin(w.*t);

figure(1)
subplot(1,2,1)
plot(t,x,'-r',t,x_exact,'ok--')
title("posição em função do tempo")
xlabel("tempo(s)")
ylabel("posição(m)")
legend('metodo de euler-implicito','solucao analitica')

subplot(1,2,2)
plot(t,v,'-r',t,v_exact,'ok--')
title("velocidade em função do tempo")
xlabel("tempo(s)")
ylabel("velocidade(m/s^2)")
legend('metodo de euler-implicito','solucao analitica')

figure(2)
plot(x,v,'.-')
title("velocidade em função a posição")
xlabel("posição(m)")
ylabel("velocidade(m/s^2)")

% Energia mecânica total (Ec+Ep)
% Ec = 1/2*m*v^2
% Ep = 0.5*K*x^2

Ec_exact = 0.5*m*(v_exact).^2;
Ep_exact = 0.5*K*(x_exact).^2;
Em_exact = Ec_exact + Ep_exact;

Ec = 0.5*m*v.^2;
Ep = 0.5*K*x.^2;
Em = Ec + Ep;

figure(3)
subplot(1,3,1)
plot(t,Ec_exact,'k.-',t,Ec,'r.-')
title("Energia Cinética")
xlabel("tempo(s)")
ylabel("Ec")
legend('solucao analitica','metodo de euler-implicito')
subplot(1,3,2)
plot(t,Ep_exact,'k.-',t,Ep,'r.-')
title("Energia potencial")
xlabel("tempo(s)")
ylabel("Ep")
legend('solucao analitica','metodo de euler-implicito')
subplot(1,3,3)
plot(t,Em_exact,'k.-',t,Em,'r.-')
title("Energia mecânica")
xlabel("tempo(s)")
ylabel("Em")
legend('solucao analitica','metodo de euler-implicito')