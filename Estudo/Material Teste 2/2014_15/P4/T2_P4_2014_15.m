%% 2o Teste Prático Física Computacional — 2014/2015 — Turma P4
%% 1 a)
clear all 
close all
clc

x0 = 0;
xf = 5;
y0 = 0;
%dy/dx = v
v0 = 0.2;

%ODE45
reltol = 1*10^(-10);
abstol_1 = 1*10^(-10);
abstol_2 = 1*10^(-10);
options = odeset('RelTol',reltol,'AbsTol',[abstol_1 abstol_2]);


[x,sol] = ode45(@f,[x0 xf],[y0 v0],options);

plot(x,sol(:,1)), xlabel('x'), ylabel('y')

%% 1 b)
clear all
close all
clc

x0 = 0;
xf = 5;
y0 = 0;

%ODE45
reltol = 1*10^(-10);
abstol_1 = 1*10^(-10);
abstol_2 = 1*10^(-10);
options = odeset('RelTol',reltol,'AbsTol',[abstol_1 abstol_2]);


% shooting
B = -2;
guess = [-0.1 0];
tol = 1E-5;
nshots = 1000;

for j = 1:nshots
    v0 = guess(j);
    [x,sol] = ode45(@f,[x0 xf],[y0 v0],options);
    result(j) = sol(end,1); 
    diff = B - result(j);
    
    if j>= 2
        m = (result(j)-result(j-1))/(guess(j)-guess(j-1));
        guess(j+1) = guess(j)+(diff)/m;
        if abs(guess(j)-guess(j-1)) < tol
            break
        end
    end
end

disp(['Valor de dy/dx inicial que faz y(x = 5) = ',num2str(result(j)),' : ',num2str(guess(j))]);
plot(x,sol(:,1)), xlabel('x'), ylabel('y')

%% 2
clear all
close all
clc

t0 = -10;
tf = 10;
N = 2^10; % so funciona com N numeros em q N é uma potencia de 2
dt = (tf-t0)/(N-1);
t = -10:dt:10;

dw = 2*pi/(N*dt);
wmin = -N/2 * dw;
wmax = (N/2 - 1)*dw;
w = wmin:dw:wmax;

y = zeros(1,N);
for i = 1:N
    if abs(t(i)) < 1
        y(i) = 1 - abs(t(i));
    end
end

Y = dt*fftshift(fft(y));

figure(1)
plot(w,abs(Y)), xlabel('w'), ylabel('abs(Y(w))')

% figure(2)
% plot(w,Y), xlabel('w'), ylabel('Y(w)')

figure(3)
hold on
plot(w,y,'k-'), xlabel('w'), ylabel('y')

% b)
yinv = ifft(ifftshift(Y))/dt;
plot(w,yinv,'r-')
legend('y(t) original', 'y(t) - Inversa Obtida')

maxdiff = max(abs(y-yinv));

disp(['Como a maior diferença ente y e yinv é quase nula ( ',num2str(maxdiff),' ) podemos concluir que o metódo tem pouco erro'])