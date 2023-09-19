
clear all
close all
clc

w = 5;
g = 9.8;
R = 0.2;

ti = 0;
tf = 10;
h = 0.1;
t = ti:h:tf;
n = length(t);

theta = zeros(1,n);         %matriz linha -> posições (n colunas)
dtheta = zeros(1,n);         %matriz linha -> velocidade (n colunas)

%condicoes iniciais
theta(1) = pi/4;
dtheta(1) = 0;

ftheta = @(dtheta) dtheta;            %derivada da posição -> velocidade
fdtheta = @(theta) (w^2*cos(theta)-g/R)*sin(theta);       %derivada da velocidade

%metodo runge kutta de quarta ordem
for i = 1:n-1
    r1v = fdtheta(theta(i));                 %sempre igual no RK
    r1x = ftheta(dtheta(i));                 %sempre igual no RK
    r2v = fdtheta(theta(i) + r1x * h/2);     %x somam com rx
    r2x = ftheta(dtheta(i) + r1v * h/2);     %v somam com rv
    r3v = fdtheta(theta(i)+ r2x * h/2);
    r3x = ftheta(dtheta(i) + r2v * h/2);
    r4v = fdtheta(theta(i) + r3x * h);
    r4x = ftheta(dtheta(i) + r3v * h);
    theta(i+1) = theta(i) + 1/6*(r1x + 2*r2x + 2*r3x + r4x)*h;
    dtheta(i+1) = dtheta(i) + 1/6*(r1v + 2*r2v + 2*r3v + r4v)*h;
end


ind = [];       %indices em que vão estar localizados os maximos
for i = 2:n-1
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

Periodo = (time(end) - time(1)) / (length(ind)-1)
%dividimos o intervalo de tempo entre os n períodos por (n-1)   


plot(t, theta, 'y')
set(gca, 'color', 'k')
title('Método Runge-Kutta 4a ordem')
xlabel 'tempo (s)'
ylabel 'theta (rad)'


  