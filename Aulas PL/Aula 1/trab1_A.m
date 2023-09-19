%% Trabalho 1 - Movimento a uma dimensão
%% A - Queda de uma pedra 
clc
clear all
close all

% f= ma <=> -mg = ma <=> (d^2/dt^2)(z) = -g     a = (d^2/dt^2)(z) uma
% dimensão(z)

% a) (d^2/dt^2)(z) = -g 
% (integrando e resolvendo vem: z = -(g/2)*t^2+v0*t+z0)

m = 150*10^-3; % massa
g = 9.8; % gravidade
h = 3; % altura de um andar
v0 = 0; % a pedra é largada
z0 = 2*h; % altura incial (segundo andar)


t = 0:0.2:1.5;
v = -g*t+v0;
z = -(g/2)*t.^2+v0*t+z0;

figure(1)
plot(t,v,'-b')
title("velocidade em função do tempo")
xlabel("tempo(s)")
ylabel("velocidade(m/s^2)")

figure(2)
plot(t,z,'-r')
title("posição em função do tempo")
xlabel("tempo(s)")
ylabel("altura(m)")

% b) Não, pois a aceleração é um equação de segunda ordem e o Metódo
% resolve equações de primeira ordem. Logo, não a podemos usar diretamente.

% c) Transforma-mos um equação linear de segunda ordem num sistema de duas
% equações lineares de primeira ordem
% dv/dt = -g e dz/dt = v (Importante!!!!!) o que leva a:
% v(k+1) = v(k) +(-g)*h <- Metódo de Euler
% z(k+1) = z(k) + v(k)*h <- Metódo de Euler

%% d)
clc
clear all
close all

g = 9.8; % gravidade
v0 = 0; % a pedra é largada
z0 = 6; % altura incial (segundo andar)

h = 0.2;
tf = 1.5;
t = 0:h:tf;

N = length(t);
v_euler = zeros(1,N);
v_euler(1) = v0;

for k = 1:(N-1) %N-1 para os comprimentos de v e t serem iguais
    v_euler(k+1) = v_euler(k)-g*h; %-g*h pois a derivada da velocidade é -g
end

% figure(1)
% plot(t,v_euler,'-b')
% title("velocidade em função do tempo")
% xlabel("tempo(s)")
% ylabel("velocidade(m/s^2)")

% e) solução analitica é a q se encontra em a)

vz = -g.*t+v0; % solução analitica

figure(2)
plot(t,v_euler,'.r-',t,vz,'ko--')
title("velocidade em função do tempo")
xlabel("tempo(s)")
ylabel("velocidade(m/s^2)")

%% f) Parte 1
clc
clear all
close all

g = 9.8; % gravidade
v0 = 0; % a pedra é largada
z0 = 6; % altura incial (segundo andar)

h = 0.2;
tf = 1.5;
t = 0:h:tf;

N = length(t);
v_euler = zeros(1,N);
v_euler(1) = v0;
z_euler = zeros(1,N);
z_euler(1) = z0;
for k = 1:(N-1) %N-1 para os comprimentos de v e t serem iguais
    v_euler(k+1) = v_euler(k)+(-g)*h; %-g*h pois a derivada da velocidade é -g
    z_euler(k+1) = z_euler(k)+v_euler(k)*h; %v(k)*h pois a derivada de z é v
end

% figure(3)
% plot(t,z_euler,'-r')
% title("posição em função do tempo")
% xlabel("tempo(s)")
% ylabel("altura(m)")

% g) solução analitica é a q se encontra em a)

z_exact = -(g/2)*t.^2+v0*t+z0; % solução analitica

figure(2)
plot(t,z_euler,'.r-',t,z_exact,'ko--')
title("posição em função do tempo")
xlabel("tempo(s)")
ylabel("altura(m)")
%% f) Parte 2 (Estimativa numérica do instante em que a pedra cai no chão e da sua velocidade.)

clc
clear all
close all

g = 9.8; % gravidade
v0 = 0; % a pedra é largada
z0 = 6; % altura incial (segundo andar)

h = 0.2;
tf = 1.5;
t = 0:h:tf;

N = length(t);
v = zeros(1,N);
v(1) = v0;
z = zeros(1,N);
z(1) = z0;
for k = 1:(N-1) %N-1 para os comprimentos de v e t serem iguais
    v(k+1) = v(k)+(-g)*h; %-g*h pois a derivada da velocidade é -g
    z(k+1) = z(k)+v(k)*h; %v(k)*h pois a derivada de z é v
    if z(k) < 0
        break %penultimo valor de z positivo e ultimo negativo, isto faz-se para se poder usar interp1
    end
end

t_impacto = interp1(z(end:-1:end-1),t(end:-1:end-1),0)   

