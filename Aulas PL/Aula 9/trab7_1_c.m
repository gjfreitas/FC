%% Trabalho Prático 7 - Transformada de Fourier discreta e sua aplicação na resolução de equações diferenciais
%% Problema 7.1: Propriedades da transformada de Fourier discreta
%% c) Depois
clear all
close all
clc

dt = 0.1;
N = 2^10*2; % 1024
t = 0:dt:(N-1)*dt;


% w com shift
% definir o eixo no espaço das frequências angulares, >=0
dw = 2*pi/(N*dt); % = delta omega
wmax = (N/2-1)*dw; % valor máximo no espaço das frequências
wmin = (-N/2)*dw;        % valor minimo no espaço das frequências
w(:,1) = wmin:dw:wmax;


% Função f(t)
y = sin(10*t) + sin(10.05*t);

% Transformada F(w)
Y1 = fft(y);
Y2 = fftshift(Y1);

% Densidade Espectral
DS2 = (dt*abs(Y2)).^2;

% Plots
figure(1)
plot(w,DS2,'.')