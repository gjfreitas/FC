%% Trabalho 1- Avançado(2020/21)
%% 1.1
% a)
clear all
close all
clc

m = 1; %kg
vlim = 6.8; % m/s
g = 9.8; % m/s^2
A = (m*g)/(vlim^2); 

h = 0.1;
t0 = 0; %s
tf = 2; %s
t = t0:h:tf; %s
N = length(t);
v = zeros(1,N);

v(1) = 0; % o volante é largado

% Metódo de Euler
% As forças que atuam sobre um volante de badminton são a força gravítica e a força de arrasto.
% Pelo enunciado, Farrasto = -A*abs(v(i))*v(i) e Fg = mg.
% Pela segunda lei de newton, assumido o eixo 'para cima' vem:
% Farrasto - Fg = ma <=> -A*abs(v(i))*v(i) - mg = ma
% Sabemos que a = (dv)/dt e divindindo a equação acima por m vem:
% (dv)/dt = (-A/m)*abs(v(i))*v(i) - g
% Então, pelo Metódo de Euler, concluimos que v(k+1) = v(k) + (-g -(A/m)*abs(v(k))*v(k))*h;

for i = 1:N-1
    v(i+1) = v(i) + (-g - (A/m)*abs(v(i))*v(i))*h;
end

% Solução Analítica
vz = -vlim*tanh((g/vlim)*t);

plot(t,vz,'-k',t,v,'m-');
xlabel("Tempo(s)")
ylabel("Velocidade(m/s)")
legend("Analitico","Euler")
str = {'h = 0.1'};
text(1.5,-3,str)
    
%% Com diferents valores de h

clear all
close all
clc

m = 1; %kg
vlim = 6.8; % m/s
g = 9.8; % m/s^2
A = (m*g)/(vlim^2); 

incr = 0.00005;

for h = 10E-4:incr:10E-1
    t = 0:h:2;
    N = length(t);
    v = zeros(1,N);
    vz = zeros(1,N);
    

    v(1) = 0; %m/s
    
    % Método de Euler
    % Metódo de Euler
    % As forças que atuam sobre um volante de badminton são a força gravítica e a força de arrasto.
    % Pelo enunciado, Farrasto = -A*abs(v(i))*v(i) e Fg = mg.
    % Pela segunda lei de newton, assumido o eixo 'para cima' vem:
    % Farrasto - Fg = ma <=> -A*abs(v(i))*v(i) - mg = ma
    % Sabemos que a = (dv)/dt e divindindo a equação acima por m vem:
    % (dv)/dt = (-A/m)*abs(v(i))*v(i) - g
    % Então, pelo Metódo de Euler, concluimos que v(k+1) = v(k) + (-g -(A/m)*abs(v(k))*v(k))*h;

    for i = 1:N-1
        v(i+1) = v(i) + (-g - (A/m)*abs(v(i))*v(i))*h;
    end
    
    % Solução Analítica
    vz = -vlim*tanh((g/vlim)*t);

    plot(t,vz,'-k',t,v,'m-');
    xlabel("Tempo(s)")
    ylabel("Velocidade(m/s)")
    legend("Analitico","Euler")
end 
    