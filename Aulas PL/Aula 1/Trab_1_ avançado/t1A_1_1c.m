%% Trabalho 1- Avançado(2020/21)
%% 1.1
% c)
clear all
close all
clc

m = 1; %kg
vlim = 6.8; % m/s
g = 9.8; % m/s^2
A = (m*g)/(vlim^2);

h = 0.01;
t = 0:h:5;
N = length(t);
v = zeros(1,N);
z = zeros(1,N);

% Condições iniciais
v0 = 16; %m/s
z0 = 1; %m

v(1) = v0;
z(1) = z0;

% Metódo de Euler
% As forças que atuam sobre um volante de badminton são a força gravítica e a força de arrasto.
% Pelo enunciado, Farrasto = -A*abs(v(i))*v(i) e Fg = mg.
% Pela segunda lei de newton, assumido o eixo 'para cima' vem:
% Farrasto - Fg = ma <=> -A*abs(v(i))*v(i) - mg = ma
% Sabemos que a = (dv)/dt e divindindo a equação acima por m vem:
% (dv)/dt = (-A/m)*abs(v(i))*v(i) - g
% Então, pelo Metódo de Euler, concluimos que v(k+1) = v(k) + (-g -(A/m)*abs(v(k))*v(k))*h;
% Sabemos também que (dz/dt) = v. Logo,pelo Metódo de Euler, z(k+1) = z(k)+v(k)*h

for i = 1:(N-1) %N-1 para os comprimentos de v e t serem iguais
    v(i+1) = v(i) + (-g - (A/m)*abs(v(i))*v(i))*h;
    z(i+1) = z(i)+v(i)*h; 
    if z(i) < 0
        break %penultimo valor de z positivo e ultimo negativo, isto faz-se para se poder usar interp1
    end
end

t_z0 = interp1(z(1:i+1),t(1:i+1),0) % Instante em que chega ao solo (s)
v_z0 = interp1(z(1:i+1),v(1:i+1),0) % Velocidade com que chega ao solo