%% 2o Teste Prático de Avaliação Discreta Física Computacional — 2015/2016 28 de abril de 2016 — Teste B
%% 1 a)
clear all
close all
clc

L = 3;
x0 = 0;
xf = L;
T = 5.0E4;
w = 1.0E5;
alfa = 5.0E-8;

h = 0.01;
x = x0:h:xf;
N = length(x);

y0 = 0;
y = nan(1,N); y(1) = y0;

%ODE45
reltol = 1*10^(-10);
abstol_1 = 1*10^(-10);
abstol_2 = 1*10^(-10);
options = odeset('RelTol',reltol,'AbsTol',[abstol_1 abstol_2]);

% shooting
B = 0;
guess = [-0.5 -0.1];
tol = 1E-12;
nshots = 100;

for j = 1:nshots
    dy0 = guess(j);
    [x,sol] = ode45(@f,[x0 xf],[y0 dy0],options,w);
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

dy_i = guess(j);
y = sol(:,1);
dy = sol(:,2);

figure(1)
plot(x,y), xlabel('x'), ylabel('y(x)'), title('Deflexão')
figure(2)
plot(x,dy), xlabel('x'), ylabel('dy'), title('Derivada da Deflexão')

disp(['dy/dx inicial => y(L) = ',num2str(result(j)),' : ',num2str(dy_i)])

%% 1 b)
clear all
close all
clc

L = 3;
x0 = 0;
xf = L;
T = 5.0E4;
alfa = 5.0E-8;

h = 0.01;
x = x0:h:xf;
N = length(x);
y0 = 0;

%ODE45
reltol = 1*10^(-10);
abstol_1 = 1*10^(-10);
abstol_2 = 1*10^(-10);
options = odeset('RelTol',reltol,'AbsTol',[abstol_1 abstol_2]);

% shooting
B = 0;
guess = [-0.5 -0.1];
tol = 1E-12;
nshots = 100;

w0 = 1E5;
wf = 5E5;
w_k = w0:100:wf;
vec = nan(length(w_k),3);
vec(:,1) = w_k(1,:);
% vec é um vetor em que a primeira coluna são os valores de w
% a segunda coluna são os valores máximos de y correspondentes a cada w
% a terceira coluna são os valores máximos de dy correspondentes a cada w


for i = 1:length(w_k)
    w = w_k(i);
    for j = 1:nshots
        dy0 = guess(j);
        [x,sol] = ode45(@f,[x0 xf],[y0 dy0],options,w);
        result(j) = sol(end,1); 
        diff = B - result(j);
        
        if j>= 2
            m = (result(j)-result(j-1))/(guess(j)-guess(j-1));
            guess(j+1) = guess(j)+(diff)/m;
            if abs(guess(j)-guess(j-1)) < tol
                    vec(i,2) = max(sol(:,1));
                    vec(i,3) = max(sol(:,2));
                break
            end
        end
    end
end

figure(1)
plot(vec(:,1),vec(:,2)), xlabel('w'), ylabel('max de y(x)'), title('Maximos da Deflexão em função de w')
figure(2)
plot(vec(:,1),vec(:,3)), xlabel('w'), ylabel('max de dy'), title('Maximos da Derivada da Deflexão em função de w')

%% 2
clear all
close all
clc

y = csvread('som2.csv');
Fs = 22.05E3; %Hz            % Frequência de Amostragem
T = 3/4;             % Periodo de Amostragem  
N = length(y);     % Length of signal
t = (0:N-1)*T;        % Time vector
dt = (t(end)-t(1))/(N-1);

f = Fs*(0:(N/2))/N;

f_sim = nan(N,1);
f_sim(1:(N/2)+1) = -f(end:-1:1);
f_sim((N/2):N) = f(1:1:end);

zf = fft(y);
z = fftshift(y);
figure(1)
plot(t,dt*abs(z)),xlabel('x'), ylabel('|Y|*dt')

dens_esp = (T*abs(z)).^2;

absf = f_sim((N/2):N);

for i = 1:length(absf)
    if absf(i) <= 1500
        freq(i) = absf(i);
        DE(i) = dens_esp(i);
    end
end

disp(['Freq 6º Harmónico: ',num2str(freq(6)),' Hz'])
figure(2)
plot(freq,DE), xlabel('freq <= 1500 HZ'), ylabel('Densidade Espectral para freq <= 1500 Hz')
