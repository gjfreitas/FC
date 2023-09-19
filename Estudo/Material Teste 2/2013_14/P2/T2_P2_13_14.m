%% 2o Teste Prático Física Computacional — 2013/2014 — Turma P2
%% 1 a)
clear all
close all
clc

x0 = -10;
xf = 10;
y0 = 10^-4;
lambda = 0.5;
v0 = y0*lambda;

%ODE45
reltol = 1*10^(-10);
abstol_1 = 1*10^(-10);
abstol_2 = 1*10^(-10);
options = odeset('RelTol',reltol,'AbsTol',[abstol_1 abstol_2]);

[x,sol] = ode45(@f,[x0 xf],[y0 v0],options,lambda);

Y = sol(:,1);
plot(x,Y), xlabel('x'), ylabel('y'), title('Solução da EDO encontrada através da rotina ODE45')

%% 1b)
clear all
close all
clc

N = 2^9;
% x = linspace(-10,10,N);
x0 = -10;
xf = 10;
h = 0.01;
x = x0:h:xf;
y0 = 10^-4;

%ODE45
reltol = 1*10^(-10);
abstol_1 = 1*10^(-10);
abstol_2 = 1*10^(-10);
options = odeset('RelTol',reltol,'AbsTol',[abstol_1 abstol_2]);

% shooting
B = 10^-4;
guess = [1.2 1.25];
tol = 1E-5;
nshots = 100;

for j = 1:nshots
    v0 = guess(j)*y0;
    [x,sol] = ode45(@f,[x0 xf],[y0 v0],options,guess(j));
    result(j) = sol(end,1); 
    diff = B - result(j);
   
    if j>= 2
        m = (result(j)-result(j-1))/(guess(j)-guess(j-1));
        guess(j+1) = guess(j)+(diff)/m;
        if abs(guess(j+1)-guess(j)) < tol
            break
        end
    end

end

lambda = guess(j);
plot(x,sol(:,1)),xlabel('x'),ylabel('y(x)'), title(strcat('lambda obtido = ',num2str(lambda)))
display(['Foi obtido um lambda = ',num2str(lambda),', para y(10) ~ ',num2str(B),', com tolerância de : ',num2str(tol)])

%% 1 d)
clear all
close all
clc

% tive que acrescentar estas linhas e mudar umas cenas para funcionar
N = 2^11;
x = linspace(-10,10,N);


y0 = 10^-4;

%ODE45
reltol = 1*10^(-10);
abstol_1 = 1*10^(-10);
abstol_2 = 1*10^(-10);
options = odeset('RelTol',reltol,'AbsTol',[abstol_1 abstol_2]);

% shooting
B = 10^-4;
guess = [1.2 1.25];
tol = 1E-5;
nshots = 100;

for j = 1:nshots
    v0 = guess(j)*y0;
    [x,sol] = ode45(@f,x,[y0 v0],options,guess(j));
    result(j) = sol(end,1); 
    diff = B - result(j);
    
    if j>= 2
        m = (result(j)-result(j-1))/(guess(j)-guess(j-1));
        guess(j+1) = guess(j)+(diff)/m;
        if abs(guess(j+1)-guess(j)) < tol
            break
        end
    end
end

dx = x(2)-x(1);
dw = 2*pi/dx;
w = -N/2*dx:dx:(N/2-1)*dx;


% Transformada de Fourier de y(x):
Y = fftshift(fft(sol(:,1)));

wmax = N*dw/2; dt = pi/wmax;
DS = (dt*abs(Y)).^2;

plot(w,DS),xlabel('x'),ylabel('(dt |y(x)|^2'), title(strcat('Densidade Espectral da solução para  ',num2str(N),' pontos'))
