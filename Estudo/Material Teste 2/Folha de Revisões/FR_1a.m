%% Folha de Revisões
clear all
close all
clc

%Constantes e Condições iniciais
m = 1.5; %kg
K = 2; %N/m
x0 = 1.9; %m
v0 = 0; %m/s
t0 = 0; %s
tf = 10; %s

alfa = -0.2;

%ODE45
reltol = 1*10^(-10);
abstol_1 = 1*10^(-10);
abstol_2 = 1*10^(-10);

options = odeset('RelTol',reltol,'AbsTol',[abstol_1 abstol_2]);
[t,sol] = ode45(@f,[t0 tf],[x0 v0],options,alfa);

x = sol(:,1);
v = sol(:,2);
    
plot(t,x), xlabel('t(s)'), ylabel('x(m)')

for k = 1:length(t)
        if x(k+1) > x(k)
            amplitude = (x(k) + x(k+1))/2;
            break
        end
end

disp(['Amplitude negativa: ',num2str(amplitude),' m'])