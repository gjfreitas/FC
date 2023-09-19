%% Trabalho 1- Avançado(2020/21)
%% 1.1
% b)
clear all
close all
clc

m = 1; %kg
vlim = 6.8; % m/s
g = 9.8; % m/s^2
A = (m*g)/(vlim^2); 

N = 6:5001;
log_error = zeros(1,length(N));
LOG_ERROR = zeros(1,length(N));
LOG_H = zeros(1,length(N));

ERROR = zeros(1,length(N));
H = zeros(1,length(N));

t = zeros(1,5001);%numero max de pontos são 5001
v = zeros(1,length(t));
vz = zeros(1,length(t));

% Condições iniciais
v(1) = 0; %m/s
    
for k = 1:length(N)
    h = 0.5/(N(k)-1);
    t(1:N(k)) = 0:h:0.5;
    
    % Metódo de Euler
    % As forças que atuam sobre um volante de badminton são a força gravítica e a força de arrasto.
    % Pelo enunciado, Farrasto = -A*abs(v(i))*v(i) e Fg = mg.
    % Pela segunda lei de newton, assumido o eixo 'para cima' vem:
    % Farrasto - Fg = ma <=> -A*abs(v(i))*v(i) - mg = ma
    % Sabemos que a = (dv)/dt e divindindo a equação acima por m vem:
    % (dv)/dt = (-A/m)*abs(v(i))*v(i) - g
    % Então, pelo Metódo de Euler, concluimos que v(k+1) = v(k) + (-g -(A/m)*abs(v(k))*v(k))*h;
    
    for i = 1:N(k)-1
        v(i+1) = v(i) + (-g - (A/m)*abs(v(i))*v(i))*h;
    end

    % Solução Analítica
    vz = -vlim*tanh((g/vlim)*t);
    
    % Cálculo de logs e erros
    error = abs(v(N(k))-vz(N(k)));
    log_error = log(error);
    log_h = log(h);
    
    LOG_ERROR(k) = log_error;
    LOG_H(k) = log_h;
    
    ERROR(k) = error;
    H(k) = h;
end

figure(1)
plot(LOG_H, LOG_ERROR,'o')
lsline;
%figure(2)
p = polyfit(LOG_H,LOG_ERROR,1);
%H1 = linspace(0,1)
%y1 = polyval(p,H1);
%plot(H,ERROR)
%hold on
%plot(H1,y1)