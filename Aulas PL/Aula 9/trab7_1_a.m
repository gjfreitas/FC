%% Trabalho Prático 7 - Transformada de Fourier discreta e sua aplicação na resolução de equações diferenciais
%% Problema 7.1: Propriedades da transformada de Fourier discreta
clear all
close all
clc

dt = 0.1;
N = 2^10; % 1024
t = 0:dt:(N-1)*dt;

% w Sem shift
dw1 = 2*pi/(N*dt); % = delta omega
wmax1 = (N-1)*dw1; % valor máximo no espaço das frequências
wmin1 = 0;        % valor minimo no espaço das frequências
w1(:,1) = wmin1:dw1:wmax1;


% w com shift
% definir o eixo no espaço das frequências angulares, >=0
dw = 2*pi/(N*dt); % = delta omega
wmax = (N/2-1)*dw; % valor máximo no espaço das frequências
wmin = (-N/2)*dw;        % valor minimo no espaço das frequências
w(:,1) = wmin:dw:wmax;


% Função f(t)
y1 = sin(t);
y2 = sin(10*t);
y3 = sin(t) + sin(10*t);

% Transformadas F(w)
Y1 = fft(y1);
Y12 = fftshift(Y1);

Y2 = fft(y2);
Y22 = fftshift(Y2);

Y3 = fft(y3);
Y32 = fftshift(Y3);


% Densidades Espectrais
DS1 = (dt*abs(Y1)).^2;      % sem shift
DS12 = (dt*abs(Y12)).^2;    % com shift

DS2 = (dt*abs(Y2)).^2;      % sem shift
DS22 = (dt*abs(Y22)).^2;    % com shift

DS3 = (dt*abs(Y3)).^2;      % sem shift
DS32 = (dt*abs(Y32)).^2;    % com shift

nyquist = pi/dt; % verificacao do numero de nyquist, que deve ser superior à menor frequência de sinal

% Plots
figure(1)
subplot(2,1,1)
plot(w1,DS1,'r.'), title('y = sin(t), sem shift')
axis([0 70 0 2000])
subplot(2,1,2)
plot(w,DS12,'r.'), title('y = sin(t), com shift')
axis([-40 40 0 2000])

figure(2)
subplot(2,1,1)
plot(w1,DS2,'b.'), title('y = sin(10*t), sem shift')
axis([0 70 0 3000])
subplot(2,1,2)
plot(w,DS22,'b.'), title('y = sin(10*t), com shift')
axis([-40 40 0 3000])

figure(3)
subplot(2,1,1)
plot(w1,DS3,'m.'), title('y = sin(t) + sin(10*t), sem shift')
axis([0 70 0 2000])
subplot(2,1,2)
plot(w,DS32,'m.'), title('y = sin(t) + sin(10*t), com shift')
axis([-40 40 0 2000])

% figure(1)
% subplot(4,2,1)
% plot(w1,DS1,'r.'), title('y = sin(t), sem shift')
% axis([0 70 0 2000])
% subplot(4,2,2)
% plot(w,DS12,'r.'), title('y = sin(t), com shift')
% axis([-40 40 0 2000])
% subplot(4,2,3)
% plot(w1,DS2,'b.'), title('y = sin(10*t), sem shift')
% axis([0 70 0 3000])
% subplot(4,2,4)
% plot(w,DS22,'b.'), title('y = sin(10*t), com shift')
% axis([-40 40 0 3000])
% subplot(4,2,5)
% plot(w1,DS3,'m.'), title('y = sin(t) + sin(10*t), sem shift')
% axis([0 70 0 2000])
% subplot(4,2,6)
% plot(w,DS32,'m.'), title('y = sin(t) + sin(10*t), com shift')
% axis([-40 40 0 2000])

