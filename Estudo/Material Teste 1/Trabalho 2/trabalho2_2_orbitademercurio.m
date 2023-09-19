%% Trabalho Prático 2
%% Problema 2.2 Órbita de Mercúrio

clear all
close all
clc

h = 0.0001; %[ano]

t(:,1) = 0:h:1;
N = length(t);

x = nan(N,1);
y = nan(N,1);
r = nan(N,1);
ang = nan(N,1);
vx = nan(N,1);
vy = nan(N,1);

x(1) = 0.47; % Unidades Astronomicas (AU)
y(1) = 0; % AU
r(1) = norm([x(1),y(1)]);
ang(1) = 0;
vx(1) = 0; % AU/ano
vy(1) = 8.2; % AU/ano

% Metódo de Euler-Cromer

for k = 1:(N-1)
    vx(k+1) = vx(k) - 4*pi^2* x(k) / (r(k)^3)*h;
    vy(k+1) = vy(k) - 4*pi^2* y(k) / (r(k)^3)*h;
    x(k+1) = x(k) + vx(k+1) * h;
    y(k+1) = y(k) + vy(k+1) * h;
    
    r(k+1) = norm([x(k+1),y(k+1)]);
    ang(k+1) = mod(atan2(y(k+1),x(k+1)),2*pi);
end

% a)
figure(1)
plot(x,y,'k.-'), xlabel('x'),ylabel('y')
axis([-0.5 0.5 -0.5 0.5])
set(gca,'PlotBoxAspectRatio',[1 1 1])

%b)

for k = 1 : N-1
    if ang(k+1) < ang(k)
        break
    end
end

disp(['periodo: ',num2str(t(k)),' anos'])
% valor observado = 88 dias ~ 0.2411 anos

% c)

index_t_meio_periodo = floor(k/2);
area = nan(index_t_meio_periodo-1,1);

for k = 1:(index_t_meio_periodo-1)
    area(k) = ( ((r(k+1) + r(k)))^2)/2*(ang(k+1)-ang(k))/2;
end

figure(2)
plot(t(1:(index_t_meio_periodo-1)),area)
title('alinea c')
xlabel('t')
ylabel('Variação da Aréa')

