%% 2o Teste Prático Física Computacional — 2013/2014 — Turma P1
%% 1 a)
clear all
close all
clc

lambda = 1;
y0 = 1;
v0 = 0;
x0 = 0;
xf = pi;

%ODE45
reltol = 1*10^(-10);
abstol_1 = 1*10^(-10);
abstol_2 = 1*10^(-10);
options = odeset('RelTol',reltol,'AbsTol',[abstol_1 abstol_2]);

[x,sol] = ode45(@f,[x0 xf],[y0 v0],options,lambda);

y = sol(:,1);
plot(x,y), xlabel('x'), ylabel('y')

%% 1 b)
clear all
close all
clc

q = 5;
y0 = 1;
v0 = 0;
x0 = 0;
xf = pi;
h = 0.01;
x = x0:h:xf;
N = length(x);

v = nan(1,N); v(1) = 0;
y = nan(1,N); y(1) = 1;

% Shooting
B = 0; % v(pi) = 0
guess = [1 1.2];

for k = 1:N-1
    v(k+1) = v(k) - ((guess(1)-q*cos(x(k)))*y(k))*h;
    y(k+1) = y(k) + v(k+1)* h ;
end
result(2) = v(end); % aqui é v pq o nosso objetivo (B) é correspondente a v

% Ja temos 2 guesses e 1 result podemos passar para o ciclo
while abs(guess(2)-guess(1)) > 1E-5 % O processo é repetido até a diferença dos guess não difira mais que uma determinada tolerância pré-estabelecida
    for k = 1:N-1
        v(k+1) = v(k) - ((guess(2)-q*cos(x(k)))*y(k))*h;
        y(k+1) = y(k) + v(k+1)* h ;
    end
    result(2) = v(end); % aqui é v pq o nosso objetivo (B) é correspondente a v
    m = (result(2)-result(1))/(guess(2)-guess(1));
    guess = [guess(2), guess(2)+ (B-result(2))/m];
    result(1) = result(2);
end

sol = guess(end-1);
disp(['lambda: ', num2str(sol)])
disp(['dy/dx (0) = ',num2str(round(result(end-1)))])
plot(x,v), xlabel('x'), ylabel('dy/dx')

%% 1 d)
clear all
close all
clc

q = 5;
y0 = 1;
v0 = 0;
x0 = 0;
xf = pi;
h = 0.01;
x = x0:h:xf;
N = length(x);

v = nan(1,N); v(1) = 0;
y = nan(1,N); y(1) = 1;

% Shooting
B = 0; % v(pi) = 0
guess = [10 10.7];

for k = 1:N-1
    v(k+1) = v(k) - ((guess(1)-q*cos(x(k)))*y(k))*h;
    y(k+1) = y(k) + v(k+1)* h ;
end
result(2) = v(end); % aqui é v pq o nosso objetivo (B) é correspondente a v

% Ja temos 2 guesses e 1 result podemos passar para o ciclo
while abs(guess(2)-guess(1)) > 1E-5 % O processo é repetido até a diferença dos guess não difira mais que uma determinada tolerância pré-estabelecida
    for k = 1:N-1
        v(k+1) = v(k) - ((guess(2)-q*cos(x(k)))*y(k))*h;
        y(k+1) = y(k) + v(k+1)* h ;
    end
    result(2) = v(end); % aqui é v pq o nosso objetivo (B) é correspondente a v
    m = (result(2)-result(1))/(guess(2)-guess(1));
    guess = [guess(2), guess(2)+ (B-result(2))/m];
    result(1) = result(2);
end
sol = guess(end-1);
disp(['lambda: ', num2str(sol)])
disp(['dy/dx (0) = ',num2str(round(result(end-1)))])

Y = fftshift(fft(y));

dx = x(2)-x(1);
dw = 2*pi/(N*dx);
w = -N/2*dw:dw:(N/2-1)*dw;

wmax = N*dw/2; dt = pi/wmax;
DS = (dt*abs(Y)).^2;
plot(w,DS), xlabel('x'), ylabel('Densidade Espectral')

