%% 1o Teste Prático Física Computacional — 2013/2014 11 março de 2014
%Turma P2

% a)
clear all
close all
clc

g = 9.8;
L = 1;
b = 1;
theta0 = 0.2;
dtheta0 = 0;
h = 0.01;
tf = 30;
t = 0:h:tf;

N = length(t);
dtheta = zeros(N,1);
theta = zeros(N,1);
theta(1) = theta0;
dtheta(1) = dtheta0;

ft = @(theta) theta;
fdt = @(theta,dtheta) -(g/L)*sin(theta)-b*dtheta;


for k = 1:N-1
    r1dt = fdt(theta(k),dtheta(k));
    r1t = ft(dtheta(k));
    
    r2dt = fdt(theta(k)+(h/2)*r1t,dtheta(k)+(h/2)*r1dt);
    r2t = ft(dtheta(k) + r1dt*(h/2));
    
    dtheta(k+1) = dtheta(k) + r2dt*h;
    theta(k+1) = theta(k) + r2t*h;
end

plot(t,theta), xlabel('t(s)'), ylabel('theta(rads)')

%% b)

clear all
close all
clc

g = 9.8;
L = 1;
b = 0.1:0.1:1;
theta0 = 0.2;
dtheta0 = 0;
h = 0.01;
tf = 30;
t = 0:h:tf;

N = length(t);
dtheta = zeros(N,1);
theta = zeros(N,1);
theta(1) = theta0;
dtheta(1) = dtheta0;

for i = 1:length(b)
    ft = @(theta) theta;
    fdt = @(theta,dtheta) -(g/L)*sin(theta)-b(i)*dtheta;


    for k = 1:N-1
        r1dt = fdt(theta(k),dtheta(k));
        r1t = ft(dtheta(k));

        r2dt = fdt(theta(k)+(h/2)*r1t,dtheta(k)+(h/2)*r1dt);
        r2t = ft(dtheta(k) + r1dt*(h/2));

        dtheta(k+1) = dtheta(k) + r2dt*h;
        theta(k+1) = theta(k) + r2t*h;
    end
    
    
    ic = 0;
    for q = 2:length(t)
        if dtheta(q-1) * dtheta(q) < 0
        ic = ic + 1;
        ind(ic,i) = q;
        end
    end
    abs_theta = abs(theta);
    for a = 1:length(ind(:,i))
        if abs_theta(ind(a,i)) < theta0/exp(1)
            tCONDICAO(a,i) = t(ind(a,i));
            break
        end
    end
end

plot(t,theta), xlabel('t(s)'), ylabel('theta(rads)')

% Instantes que correspondem à condição pedida.
% cada linha corresponde a um 'b' diferente, variando este entre 0,1 até 1,
% com um passo de 0,1.
t_resposta = nonzeros(tCONDICAO)

%% c)
clear all
close all
clc

g = 9.8;
L = 1;
b = 1;
theta0 = 0.2;
dtheta0 = 0;
h = 0.01;
tf = 30;
t = 0:h:tf;

N = length(t);
dtheta = zeros(N,1);
theta = zeros(N,1);
theta(1) = theta0;
dtheta(1) = dtheta0;
omega = 1;

ft = @(theta) theta;
fdt = @(t, theta, dtheta) -(g/L)*sin(theta)-b*dtheta + 10 * sin(omega * t);


for k = 1:N-1
    r1dt = fdt(t(k), theta(k),dtheta(k));
    r1t = ft(dtheta(k));
    
    r2dt = fdt(t(k)+ h/2, theta(k)+(h/2)*r1t,dtheta(k)+(h/2)*r1dt);
    r2t = ft(dtheta(k) + r1dt*(h/2));
    
    dtheta(k+1) = dtheta(k) + r2dt*h;
    theta(k+1) = theta(k) + r2t*h;
end
hold on
plot(t,theta), xlabel('t(s)'), ylabel('theta(rads)')

% Indexação dos pontos onde a velocidade troca de sentido e interpolação:
ic = 0;
for n = floor(length(t)/2):length(t)
    if dtheta(n-1)*dtheta(n)<0
        ic = ic+1;
        ind(ic) = n;
    end
end

for n = 1:ic
    tI1(n,:) = t(ind(n)-1:ind(n)+1);
    thetaI1(n,:) = theta(ind(n)-1:ind(n)+1);
    dthetaI1(n,:) = dtheta(ind(n)-1:ind(n)+1);
end
for n = 1:ic
    tN1(n) = interp1(dthetaI1(n,:),tI1(n,:),0);
    thetaN1(n) = interp1(dthetaI1(n,:),thetaI1(n,:),0);
end

% Cálculo do período através da média do tempo entre picos.
Periodo1 = mean(diff(tN1))*2;

% Cálculo da amplitude através da média da posição entre extremos positivos
% e negativos.
for j = 1:length(thetaN1)/2
    thetaSN1(j) = abs(thetaN1(2*j)-abs(thetaN1(2*j-1)));
end
Amplitude1 = mean(thetaSN1)/2;

disp(['Periodo para omega = 1 : ',num2str(Periodo1), ' s'])
disp(['Amplitude para omega = 1 : ',num2str(Amplitude1), ' m'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
% Omega = 2;

g = 9.8;
L = 1;
b = 1;
theta0 = 0.2;
dtheta0 = 0;
h = 0.01;
tf = 30;
t = 0:h:tf;

N = length(t);
dtheta = zeros(N,1);
theta = zeros(N,1);
theta(1) = theta0;
dtheta(1) = dtheta0;
omega = 2;

ft = @(theta) theta;
fdt = @(t, theta, dtheta) -(g/L)*sin(theta)-b*dtheta + 10 * sin(omega * t);


for k = 1:N-1
    r1dt = fdt(t(k), theta(k),dtheta(k));
    r1t = ft(dtheta(k));
    
    r2dt = fdt(t(k)+ h/2, theta(k)+(h/2)*r1t,dtheta(k)+(h/2)*r1dt);
    r2t = ft(dtheta(k) + r1dt*(h/2));
    
    dtheta(k+1) = dtheta(k) + r2dt*h;
    theta(k+1) = theta(k) + r2t*h;
end

plot(t,theta), xlabel('t(s)'), ylabel('theta(rads)')

% Indexação dos pontos onde a velocidade troca de sentido e interpolação:
ic = 0;
for n = floor(length(t)/2):length(t)
    if dtheta(n-1)*dtheta(n)<0
        ic = ic+1;
        ind(ic) = n;
    end
end
for n = 1:ic
    tI1(n,:) = t(ind(n)-1:ind(n)+1);
    thetaI1(n,:) = theta(ind(n)-1:ind(n)+1);
    dthetaI1(n,:) = dtheta(ind(n)-1:ind(n)+1);
end
for n = 1:ic
    tN1(n) = interp1(dthetaI1(n,:),tI1(n,:),0);
    thetaN1(n) = interp1(dthetaI1(n,:),thetaI1(n,:),0);
end

% Cálculo do período através da média do tempo entre picos.
Periodo2 = mean(diff(tN1))*2;

% Cálculo da amplitude através da média da posição entre extremos positivos
% e negativos.
for j = 1:length(thetaN1)/2
    thetaSN1(j) = abs(thetaN1(2*j)-abs(thetaN1(2*j-1)));
end
Amplitude2 = mean(thetaSN1)/2;

disp(['Periodo para omega = 2 : ',num2str(Periodo2), ' s'])
disp(['Amplitude para omega = 2 : ',num2str(Amplitude2), ' m'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
% Omega = 3;

g = 9.8;
L = 1;
b = 1;
theta0 = 0.2;
dtheta0 = 0;
h = 0.01;
tf = 30;
t = 0:h:tf;

N = length(t);
dtheta = zeros(N,1);
theta = zeros(N,1);
theta(1) = theta0;
dtheta(1) = dtheta0;
omega = 1;

ft = @(theta) theta;
fdt = @(t, theta, dtheta) -(g/L)*sin(theta)-b*dtheta + 10 * sin(omega * t);


for k = 1:N-1
    r1dt = fdt(t(k), theta(k),dtheta(k));
    r1t = ft(dtheta(k));
    
    r2dt = fdt(t(k)+ h/2, theta(k)+(h/2)*r1t,dtheta(k)+(h/2)*r1dt);
    r2t = ft(dtheta(k) + r1dt*(h/2));
    
    dtheta(k+1) = dtheta(k) + r2dt*h;
    theta(k+1) = theta(k) + r2t*h;
end

plot(t,theta), xlabel('t(s)'), ylabel('theta(rads)')
legend('Omega = 1','Omega = 2','Omega = 3')

% Indexação dos pontos onde a velocidade troca de sentido e interpolação:
ic = 0;
for n = floor(length(t)/2):length(t)
    if dtheta(n-1)*dtheta(n)<0
        ic = ic+1;
        ind(ic) = n;
    end
end
for n = 1:ic
    tI1(n,:) = t(ind(n)-1:ind(n)+1);
    thetaI1(n,:) = theta(ind(n)-1:ind(n)+1);
    dthetaI1(n,:) = dtheta(ind(n)-1:ind(n)+1);
end
for n = 1:ic
    tN1(n) = interp1(dthetaI1(n,:),tI1(n,:),0);
    thetaN1(n) = interp1(dthetaI1(n,:),thetaI1(n,:),0);
end

% Cálculo do período através da média do tempo entre picos.
Periodo3 = mean(diff(tN1))*2;

% Cálculo da amplitude através da média da posição entre extremos positivos
% e negativos.
for j = 1:length(thetaN1)/2
    thetaSN1(j) = abs(thetaN1(2*j)-abs(thetaN1(2*j-1)));
end
Amplitude3 = mean(thetaSN1)/2;

disp(['Periodo para omega = 3 : ',num2str(Periodo3), ' s'])
disp(['Amplitude para omega = 3 : ',num2str(Amplitude3), ' m'])

