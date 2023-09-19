%% Exame Prático — Física Computacional — 2018/2019 21 de junho de 2019
% Apenas o ex 1 sai para o teste prático 1 2020/21

% 1 a)
clear all
close all
clc

vlim = 6.8; %m/s
z0 = 1; %m
x0 = 0; %m
v0 = 30;%m/s
v0x = v0 * cosd(60);
v0z = v0 * sind(60);
g = 9.8;
t0 = 0;
tf = 10;
h = 0.01;
t = t0:h:tf;

N = length(t);
vx = nan(1,N); vx(1) = v0x;
vz = nan(1,N); vz(1) = v0z;
z = nan(1,N); z(1) = z0;
x = nan(1,N); x(1) = x0;

fz = @(vz) vz;
fvz = @(vz) (- g/(vlim)^2 * vz*abs(vz) - g);
fx = @(vx) vx;
fvx = @(vx) (- g/(vlim)^2 * vx*abs(vx));

% Runge Kutta 2a ordem
for i = 1:N-1
    r1vx = fvx(vx(i));
    r1vz = fvz(vz(i));
    r1z = fz(vz(i));
    r1x = fx(vx(i));
    
    r2vx = fvx(vx(i)+r1vx*h);
    r2vz = fvz(vz(i)+r1vz*h);
    r2z = fz(vz(i)+r1vz*h);
    r2x = fx(vx(i)+r1vx*h);
    
    vx(i+1) = vx(i) + r1vx*h/2 + r2vx*h/2;
    vz(i+1) = vz(i) + r2vx*h/2 + r2vz*h/2;
    x(i+1) = x(i) + r1x*h/2 + r2x*h/2;
    z(i+1) = z(i) + r1z*h/2 + r2z*h/2;
    
    %parar quando atinge o solo
    if z(i) < 0
        break
    end
end 

figure(1)
plot(x,z,'m'), xlabel('x (m)'), ylabel('altura (m)')
title('Trajetória do volante')

figure(2)
plot(t,z,'r'), xlabel('tempo (s)'), ylabel('altura (m)')
title('Evolução da altura do volante')


%interpolacao para descobrir o ponto em que a função z(t) é 0
a = interp1(z(1:i),t(1:i),0);
disp(['Tempo de voo : ', num2str(a), ' segundos'])


%% 1 b)
clear all
close all
clc


z0 = 1; %m
x0 = 0; %m
v0 = 30;%m/s
vx0 = v0 * cosd(60);
vz0 = v0 * sind(60);

t0 = 0;
tf = 2.5;
h = 0.01;
t = t0:h:tf;

%tolerancias dadas no enunciado
reltol = 3E-14;
abstol_1=1E-13;
abstol_2=1E-13;


options = odeset('RelTol', reltol, 'AbsTol', [abstol_1 abstol_2 abstol_1 abstol_2]);      %opcoes para o ode45
[t, sol] = ode45(@f2, t, [x0 vx0 z0 vz0], options);              %usa a funcao f em ficheiro separado
%devolve um vetor t e um vetor sol      (tempo = t)
%a primeira coluna de sol tem todas as posições      (x = sol(:,1)
%a segunda coluna de sol tem todas as velocidades     (vx = sol(:,2)

figure(1)
plot(sol(:,1), sol(:,3), '-m'), xlabel('Posição (m)'), ylabel('Altura (m)')
title('Volante de badminton - ode45')


function derivadas = f2(t,solucao)
    g = 9.8;
    vlim = 6.8;
    derivadas = zeros(4,1);
    % Esta linha é necessária e faz do vetor de saída um vetor coluna.

    % O vetor solucao tem os valores de x e v para o tempo t em que a
    % função é chamada pela rotina ode45.
    derivadas(1) = solucao(2);
    % Escreva acima a expressão da derivada de x em função de solucao(1)
    % e de solucao(2).
    derivadas(2) = (- g/(vlim)^2 * solucao(2)*abs(solucao(2)));
    % Escreva a acima expressão da derivada de v em função de solucao(1)
    % e de solucao(2).
    derivadas(3) = solucao(4);
    derivadas(4) = (- g/(vlim)^2 * solucao(4)*abs(solucao(4)) - g);
    % x = solucao(1);
    % vx = solucao(2);
    % z = solucao(3);
    % vz = solucao(4);

end
