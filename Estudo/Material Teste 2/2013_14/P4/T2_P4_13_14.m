%% 2o Teste Prático Física Computacional — 2013/2014 — Turma P4
%% 1 a)
clear all
close all
clc

x0 = 0;
xf = 5;
u0 = 1;
Q0 = -0.1;  % du/dt = Q

%ODE45
reltol = 1*10^(-10);
abstol_1 = 1*10^(-10);
abstol_2 = 1*10^(-10);
options = odeset('RelTol',reltol,'AbsTol',[abstol_1 abstol_2]);


[x,sol] = ode45(@f,[x0 xf],[u0 Q0],options);
u = sol(:,1);
plot(x,u), xlabel('x(m)'), ylabel('u')

%% 1 b)
clear all
close all
clc

const = 1;
x0 = 0;
xf = 5;
h = 0.01;
x = x0:h:xf;
N = length(x);

u = nan(1,N); u(1) = 1;

% Shooting
guess = [-1.0 -0.5];
B = 1.5;
nshots = 10000;

for ishot = 1:nshots
     Q = nan(1,N);
     Q(1) = guess(ishot);
     for k = 1:N-1
         Q(k+1) = Q(k) + (const*sqrt(1 + (Q(k))^2))*h;
         u(k+1) = u(k) + Q(k+1) *h;
     end
     result(ishot) = u(end);

     if ishot >= 2
        m = (result(ishot)-result(ishot-1))/(guess(ishot)-guess(ishot-1));
        guess(ishot+1) = guess(ishot) + (B - result(ishot))/m;
        if abs(guess(ishot) - guess(ishot-1)) < 1E-6
            break
        end   
     end
end

solucao = guess(ishot);
disp(['u(5) : ',num2str(result(ishot)),' m'])
disp(['Derivada de u(0) : ',num2str(solucao)])
plot(x,u), xlabel('x'), ylabel('u')

%% 2 provavelmente ta mal mas ya
clear all
close all
clc

% considerando N = 1024 = 2^10 e dw = 0.1;
N = 2^10;

% w com shift
% definir o eixo no espaço das frequências angulares, >=0
dw = 0.1; % = delta omega
wmax = (N/2-1)*dw; % valor máximo no espaço das frequências
wmin = (-N/2)*dw;        % valor minimo no espaço das frequências
w = wmin:dw:wmax;

dt = 2*pi/(N * dw);
t = 0:dt: (N-1)*dt;

y = cos(140*t)+cos(140.5*t);

Y = fftshift(fft(y));
DS = (dt*abs(Y)).^2;    % com shift

plot(w,DS)
