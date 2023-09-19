%% Exame Prático de Recurso — Parte 1 Física Computacional — 2018/2019 9 de julho de 2019
clear all
close all
clc

w = 5; 
% b) w = 8; w = 9; 
% c) w = 9;
g = 9.8;
R = 0.2;
theta0 = pi/4;
% c) theta0 = acos(g/(w^2 * R));
dtheta0 = 0;


ti = 0;
tf = 10;
h = 0.1;
t = ti:h:tf;
N = length(t);

theta = zeros(1,N);
dtheta = zeros(1,N);

theta(1) = theta0;
dtheta(1) = dtheta0;

ftheta = @(dtheta) dtheta;
fdtheta = @(theta) (w^2*cos(theta)-g/R)*sin(theta);

% Runge Kutta 4
for i = 1:N-1
    r1v = fdtheta(theta(i));  
    r1x = ftheta(dtheta(i));
    
    r2v = fdtheta(theta(i) + r1x * h/2); 
    r2x = ftheta(dtheta(i) + r1v * h/2);
    
    r3v = fdtheta(theta(i)+ r2x * h/2);
    r3x = ftheta(dtheta(i) + r2v * h/2);
    
    r4v = fdtheta(theta(i) + r3x * h);
    r4x = ftheta(dtheta(i) + r3v * h);
    
    theta(i+1) = theta(i) + 1/6*(r1x + 2*r2x + 2*r3x + r4x)*h;
    dtheta(i+1) = dtheta(i) + 1/6*(r1v + 2*r2v + 2*r3v + r4v)*h;
end


ind = [];       %indices em que vão estar localizados os maximos
for i = 2:N-1
    if theta(i) >= theta(i-1) && theta(i) >= theta(i+1)     %condicao de maximo
        ind = [ind i];                      %adiciona o novo indice ao fim do vetor
    end
end

max = [];                       %valores dos maximos
time = [];                      %valores do tempo em que se regista o maximo

for i = 1:length(ind)       %percorre todos os indices em que vao estar localizados maximos
    k = ind(i);                 %indice do maximo com que estamos a trabalhar
   
    aux = lagr (t(k-1:k+1), theta(k-1:k+1));       
    %usa o script lagr.m -> devolve [tempo do maximo; maximo] ajustados por lagrange
                                             
    max = [max  aux(2) ];
    time = [time aux(1) ];
end


periodo = (time(end) - time(1)) / (length(ind)-1);
%dividimos o intervalo de tempo entre os n períodos por (n-1) 
disp(['Periodo : ',num2str(periodo),' s'])


plot(t, theta)
title('Método Runge-Kutta 4a ordem')
xlabel('tempo (s)')
ylabel('theta (rad)')

%% d) Condiões de c) usando ode45
  clear all
close all
clc

w = 9;
g = 9.8;
R = 0.2;
theta0 = acos(g/(w^2 * R));
dtheta0 = 0;


ti = 0;
tf = 10;
h = 0.1;
t = ti:h:tf;

reltol = 3E-14;
abstol_1=1E-13;
abstol_2=1E-13;

options = odeset('RelTol', reltol, 'AbsTol', [abstol_1 abstol_2]);      %opcoes para o ode45
[t_ode45, sol] = ode45(@f, t, [theta0 dtheta0], options);

theta = sol(:,1);

plot(t_ode45, theta, 'y')
plot(t, theta, 'm')
title('Método Runge-Kutta 4a ordem')
xlabel('tempo (s)')
ylabel('theta (rad)')