%% Folha de Revisões
clear all
close all
clc

%Constantes e Condições iniciais
m = 1.5; %kg
K = 2; %N/m
x0 = 1.9; %m
v0 = 0; %m/s
alfa = -0.2;
t0 = 0;
h = 0.01;
tf = 10;

%ODE45
reltol = 1*10^(-10);
abstol_1 = 1*10^(-10);
abstol_2 = 1*10^(-10);
options = odeset('RelTol',reltol,'AbsTol',[abstol_1 abstol_2]);

% objetivo
B = -1.5; %m (amplitude)

%tolerancia
tol = 1E-3;

% guesses iniciais para alfa
guess(1) = -0.2;
guess(2) = -0.25;
result = [nan,nan];

[t,sol] = ode45(@f,[t0 tf],[x0 v0],options,guess(1));
x = sol(:,1);

for k = 1:length(t)
    if x(k+1) > x(k)
        result(1) = (x(k) + x(k+1))/2;
        break
    end
end


% Ja temos 2 guesses e 1 result podemos passar para o ciclo
while abs(guess(2)-guess(1)) > tol
    clear t
    clear x
    
    [t,sol] = ode45(@f,[t0 tf],[x0 v0],options,guess(2));
    
    x = sol(:,1);
    
    for k = 1:length(t)
        if x(k+1) > x(k)
            result(2) = (x(k) + x(k+1))/2;
            break
        end
    end
    
    m = (result(2)-result(1))/(guess(2)-guess(1));
    guess(3) = guess(2)+ (B-result(2))/m;
    
    guess(1) = guess(2);
    result(1) = result(2);
    
    guess(2) = guess(3);
    
end

sol = guess(end-1);
disp(['alfa: ', num2str(sol)])
figure(2)
plot(t,x,'.-'), xlabel('t(s)'),ylabel('x(m)')
grid