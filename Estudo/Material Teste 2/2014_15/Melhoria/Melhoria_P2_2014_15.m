%% Teste Prático de Melhoria Parte 2 Física Computacional — 2014/2015
%% 1 a)
clear all
close all
clc


t0 = 0;
tf = 50;
theta0 = 0.2;
% dtheta/dt = Q
Q0 = 0;


%ODE45
reltol = 1*10^(-10);
abstol_1 = 1*10^(-10);
abstol_2 = 1*10^(-10);
options = odeset('RelTol',reltol,'AbsTol',[abstol_1 abstol_2]);
[t,sol] = ode45(@f,[t0 tf],[theta0 Q0],options);
theta = sol(:,1);


% Euler-Cromer
% h = 0.01;
% t = t0:h:tf;
% N = length(t);
% theta = nan(1,N); theta(1) = theta0;
% Q = nan(1,N); Q(1) = Q0;
% 
% for k = 1:N-1
%     Q(k+1) = Q(k) + (-sin(theta(k)))*h;
%     theta(k+1) = theta(k) + Q(k+1)*h;
% end

plot(t,theta), xlabel('t'), ylabel('theta')

% b)
j = 1;
for i = 2 : length(theta)-1
    if theta(i) >= theta(i+1) && theta(i) >= theta(i-1)
        maxTheta(j) = i;
        j = j+1;
    end
end    

periodo = (t(maxTheta(2))-t(maxTheta(1))); % periodo pode ser calculado como a diferença de dois máximos

%Outra forma de calcular o periodo
[pks,locs] = findpeaks(theta);                                     % Peaks & Locations
mean_period = mean(diff(t(locs))); 

freq = 1/periodo;
disp(['Frequência: ',num2str(freq),' rad/s'])


%% 1 c)
clear all
close all
clc

%ODE45
reltol = 1*10^(-10);
abstol_1 = 1*10^(-10);
abstol_2 = 1*10^(-10);
options = odeset('RelTol',reltol,'AbsTol',[abstol_1 abstol_2]);


t0 = 0;
tf = 50;
% dtheta/dt = Q
Q0 = 0;

h = 0.01;
t = t0:h:tf;
N = length(t);
theta = nan(1,N);
Q = nan(1,N); Q(1) = Q0;

% shooting
B = 1;
guess = [0 0.2];    % guesses para theta0
nshots = 500;
tol = 1E-6;

for j = 1:nshots
    clear theta
    clear period
    
    [t,sol] = ode45(@f,[t0 tf],[guess(j) Q0],options);
    theta = sol(:,1);
    u = 1;
    for i = 2 : length(theta)-1
        if theta(i) >= theta(i+1) && theta(i) >= theta(i-1)
            maxTheta(u) = i;
            u = u+1;
        end
    end    

    period = abs((t(maxTheta(2))-t(maxTheta(1)))); 
    result(j) = 1/period;    % o nosso result é a freq
    diff = B - result(j);
    
    if j>= 2
        m = (result(j)-result(j-1))/(guess(j)-guess(j-1));
        guess(j+1) = guess(j)+(diff)/m;
        if abs(guess(j)-guess(j-1)) < tol
            break
        end
    end
end

T0_sol = guess(j);
freq = result(j);
disp(['Theta(0) => Freq = ',num2str(freq),' Hz : ',num2str(T0_sol)])

%% 1 d) provavelmente ta mal mas ya
clearvars -except T0_sol
close all
clc

theta0 = 0.092367289641677; %T0_sol
% theta0 = T0_sol;
t0 = 0;
tf = 50;
% dtheta/dt = Q
Q0 = 0;

%ODE45
reltol = 1*10^(-10);
abstol_1 = 1*10^(-10);
abstol_2 = 1*10^(-10);
options = odeset('RelTol',reltol,'AbsTol',[abstol_1 abstol_2]);
[t,sol] = ode45(@f,[t0 tf],[theta0 Q0],options);
theta = sol(:,1);

N = length(theta);
dt = (t(end)-t(1))/(N-1);
THETA = dt * fftshift(fft(theta)); % tem que se multiplicar por dt para obter a transformada de Fourier pq o matlab calcula a Transformada de Fourier Discreta

dw = 2*pi/(N*dt);
wmin = -N/2 * dw;
wmax = (N/2 - 1)*dw;
w = wmin:dw:wmax;

Y = dt*fftshift(fft(theta));
DE = (dt*abs(Y)).^2;

plot(w,DE);

for i = 2:N-1
    if abs(THETA(i)) >= abs(THETA(i+1)) && abs(THETA(i)) >= abs(THETA(i-1))
        maxT = THETA(i);
        ind = i;
    end
end

disp(['w correspondete ao máximo de THETA(w) : ',num2str(w(ind))])

%% 1 e) n sei fazer pelo metódo das transformadas pq n me lembro de ter falado disso
clear all
close all
clc

theta0_k = 0.2:0.1:1;
t0 = 0;
tf = 50;
% dtheta/dt = Q
Q0 = 0;

%ODE45
reltol = 1*10^(-10);
abstol_1 = 1*10^(-10);
abstol_2 = 1*10^(-10);
options = odeset('RelTol',reltol,'AbsTol',[abstol_1 abstol_2]);

for j = 1:length(theta0_k)
    theta0 = theta0_k(j);
    [t,sol] = ode45(@f,[t0 tf],[theta0 Q0],options);
    theta = sol(:,1);
    u = 1;
    for i = 2 : length(theta)-1
        if theta(i) >= theta(i+1) && theta(i) >= theta(i-1)
            maxTheta(u) = i;
            u = u+1;
        end
    end    

    period = abs((t(maxTheta(2))-t(maxTheta(1)))); 
    freq(j) = 1/period;
    
end

plot(theta0_k,freq), xlabel('Theta(0)'), ylabel('Frequência'), title('Frequência em função de Theta(0)')