%% Melhoria do 1o Teste Prático Física Computacional — 2015/2016 15 de junho de 2016

clear all 
close all
clc

% 1 a) Um Fd de cada vez
w0 = 1;
q = 1/2;
wd = 2/3;
theta0 = 0.2;
j0 = 0; % seja j a derivada de theta em ordem ao tempo

Fd = 0;
% Fd = 0.1;
% Fd = 1.2;
% Fd = [ 0 0.1 1.2];

h = 0.01;
t = 0:h:100;
N = length(t);
j = nan(1,N); j(1) = j0;
theta = nan(1,N); theta(1) = theta0;

fj = @(t, theta, j) -w0 * sin(theta) - q * j + Fd * sin(wd * t);
ftheta = @(t,j) j;

for k = 1:N-1
    r1theta = ftheta(t(k), j(k));
    r1j = fj(t(k), theta(k), j(k));
    
    r2theta = ftheta(t(k)+h/2, j(k)+r1j*h/2); 
    r2j = fj(t(k)+h/2, theta(k)+r1theta*h/2, j(k)+r1j*h/2);
        
    r3theta = ftheta(t(k)+h/2, j(k)+r2j*h/2); 
    r3j = fj(t(k)+h/2, theta(k)+r2theta*h/2, j(k)+r2j*h/2);
    
    r4theta = ftheta(t(k)+h, j(k)+r3j*h); 
    r4j = fj(t(k)+h, theta(k)+r3theta*h, j(k)+r3j*h);
    
    theta(k+1) = theta(k) + 1/6 *(r1theta + 2*r2theta + 2*r3theta +r4theta)*h;
    j(k+1) = j(k) + 1/6 *(r1j + 2*r2j + 2*r3j +r4j)*h;
end

plot(t,theta), xlabel('t(s)'), ylabel('theta')


%% 1 a) de outra forma -> Todos os Fd no mesmo gráfico

clear all 
close all
clc

% 1 a)
w0 = 1;
q = 1/2;
wd = 2/3;
theta0 = 0.2;
j0 = 0; % seja j a derivada de theta em ordem ao tempo

% Fd = 0;
% Fd = 0.1;
% Fd = 1.2;
Fd = [0, 0.1, 1.2];

for i = 1:length(Fd)
    h = 0.01;
    t = 0:h:100;
    N = length(t);
    j = nan(1,N); j(1) = j0;
    theta = nan(1,N); theta(1) = theta0;

    fj = @(t, theta, j) -w0 * sin(theta) - q * j + Fd(i) * sin(wd * t);
    ftheta = @(t,j) j;

    for k = 1:N-1
        r1theta = ftheta(t(k), j(k));
        r1j = fj(t(k), theta(k), j(k));

        r2theta = ftheta(t(k)+h/2, j(k)+r1j*h/2); 
        r2j = fj(t(k)+h/2, theta(k)+r1theta*h/2, j(k)+r1j*h/2);

        r3theta = ftheta(t(k)+h/2, j(k)+r2j*h/2); 
        r3j = fj(t(k)+h/2, theta(k)+r2theta*h/2, j(k)+r2j*h/2);

        r4theta = ftheta(t(k)+h, j(k)+r3j*h); 
        r4j = fj(t(k)+h, theta(k)+r3theta*h, j(k)+r3j*h);

        theta(k+1) = theta(k) + 1/6 *(r1theta + 2*r2theta + 2*r3theta +r4theta)*h;
        j(k+1) = j(k) + 1/6 *(r1j + 2*r2j + 2*r3j +r4j)*h;
    end

    figure(1)
    hold on
    
    if Fd(i) == 0
        plot(t,theta,'b-'), xlabel('t(s)'), ylabel('theta')
    end

    if Fd(i) == 0.1
        plot(t,theta,'m-'), xlabel('t(s)'), ylabel('theta')
    end
    
    if Fd(i) == 1.2
        plot(t,theta,'k-'), xlabel('t(s)'), ylabel('theta')
    end
    
   
end
legend('Fd = 0','Fd = 0.1','Fd = 1.2')

% 1 d) obtivemos uma trajetoria ....

%% 1 b) Fd = 0  -> maximos de theta
clear all
close all
clc

w0 = 1;
q = 1/2;
wd = 2/3;
theta0 = 0.2;
j0 = 0; % seja j a derivada de theta em ordem ao tempo

Fd = 0;

h = 0.01;
t = 0:h:100;
N = length(t);
j = nan(1,N); j(1) = j0;
theta = nan(1,N); theta(1) = theta0;

fj = @(t, theta, j) -w0 * sin(theta) - q * j + Fd * sin(wd * t);
ftheta = @(t,j) j;

for k = 1:N-1
    r1theta = ftheta(t(k), j(k));
    r1j = fj(t(k), theta(k), j(k));
    
    r2theta = ftheta(t(k)+h/2, j(k)+r1j*h/2); 
    r2j = fj(t(k)+h/2, theta(k)+r1theta*h/2, j(k)+r1j*h/2);
        
    r3theta = ftheta(t(k)+h/2, j(k)+r2j*h/2); 
    r3j = fj(t(k)+h/2, theta(k)+r2theta*h/2, j(k)+r2j*h/2);
    
    r4theta = ftheta(t(k)+h, j(k)+r3j*h); 
    r4j = fj(t(k)+h, theta(k)+r3theta*h, j(k)+r3j*h);
    
    theta(k+1) = theta(k) + 1/6 *(r1theta + 2*r2theta + 2*r3theta +r4theta)*h;
    j(k+1) = j(k) + 1/6 *(r1j + 2*r2j + 2*r3j +r4j)*h;
end
figure(1)
plot(t,theta), xlabel('t(s)'), ylabel('theta')


count = 0;

for k = 2:(N-1) % começamos em 2 pois se começarmos em 1 não temos nenhum anterior para comparar (em matlab não ha posição 0)
    if theta(k-1) <= theta(k) && theta(k) >= theta(k+1)
        count = count + 1;
        
        aux = lagr(t(k-1:k+1),theta(k-1:k+1));
        t_max(count) = aux(1);
        theta_max(count) = aux(2);
    end
end
figure(2)
plot(t_max,log(theta_max), '*b'), xlabel('tmax'), ylabel('log(theta_max)')
lsline
p = polyfit(t_max,log(theta_max), 1);
a = p(1); %declive
b = p(2); % ordenada na origem
disp(['Theta_max =  ',num2str(b),' - (',num2str(a),')*t_max'])
disp('a = -q/2')


%% 1 c) Fd = 0.1 frequencia de oscilação apos t = 50
clear all
close all
clc

w0 = 1;
q = 1/2;
wd = 2/3;
theta0 = 0.2;
j0 = 0; % seja j a derivada de theta em ordem ao tempo

Fd = 0.1;

h = 0.01;
t = 0:h:50;
N = length(t);
j = nan(1,N); j(1) = j0;
theta = nan(1,N); theta(1) = theta0;

fj = @(t, theta, j) -w0 * sin(theta) - q * j + Fd * sin(wd * t);
ftheta = @(t,j) j;

for k = 1:N-1
    r1theta = ftheta(t(k), j(k));
    r1j = fj(t(k), theta(k), j(k));
    
    r2theta = ftheta(t(k)+h/2, j(k)+r1j*h/2); 
    r2j = fj(t(k)+h/2, theta(k)+r1theta*h/2, j(k)+r1j*h/2);
        
    r3theta = ftheta(t(k)+h/2, j(k)+r2j*h/2); 
    r3j = fj(t(k)+h/2, theta(k)+r2theta*h/2, j(k)+r2j*h/2);
    
    r4theta = ftheta(t(k)+h, j(k)+r3j*h); 
    r4j = fj(t(k)+h, theta(k)+r3theta*h, j(k)+r3j*h);
    
    theta(k+1) = theta(k) + 1/6 *(r1theta + 2*r2theta + 2*r3theta +r4theta)*h;
    j(k+1) = j(k) + 1/6 *(r1j + 2*r2j + 2*r3j +r4j)*h;
end
figure(1)
plot(t,theta), xlabel('t(s)'), ylabel('theta')

count = 0;

for k = 2:(N-1) % começamos em 2 pois se começarmos em 1 não temos nenhum anterior para comparar (em matlab não ha posição 0)
    if theta(k-1) <= theta(k) && theta(k) >= theta(k+1)
        count = count + 1;
        
        aux = lagr(t(k-1:k+1),theta(k-1:k+1));
        t_max(count) = aux(1);

    end
end

p = mean(diff(t_max)) ; %periodo
disp(['Período médio: ',num2str(p),' s'])
freq = 1/p;
disp(['Frequência =  ',num2str(freq),' Hz '])
i = freq/wd;
disp(['Frequência =   ',num2str(i),' * Wd'])

%% 2 Euler implicito
clear all
close all
clc

w0 = 1;
q = 1/2;
wd = 2/3;
theta0 = 0.2;
j0 = 0; % seja j a derivada de theta em ordem ao tempo

Fd = 1.2;
t0 = 0;

h = 0.01;
t = t0:h:100;
N = length(t);
j = nan(1,N); j(1) = j0;
theta = nan(1,N); theta(1) = theta0;

A = [1, -h; w0 * h, 1 + q * h];
b = [theta0; j0 + Fd * h * sin(wd * t0)];
for k = 1:(N-1)
    z = linsolve(A,b);
    theta(k+1) = z(1);
    j(k+1) = z(2);
    
    b = [z(1); z(2) + Fd * h * sin(wd * t(k))];
end

figure(1)
plot(t,theta), xlabel('t(s)'), ylabel('theta')

count = 0;

for k = 2:(N-1) % começamos em 2 pois se começarmos em 1 não temos nenhum anterior para comparar (em matlab não ha posição 0)
    if theta(k-1) <= theta(k) && theta(k) >= theta(k+1)
        count = count + 1;
        
        aux = lagr(t(k-1:k+1),theta(k-1:k+1));
        t_max(count) = aux(1);
    end
end

p = mean(diff(t_max)) ; %periodo
disp(['Período médio: ',num2str(p),' s'])
freq = 1/p;
disp(['Frequência =  ',num2str(freq),' Hz '])
i = freq/wd;
disp(['Frequência =   ',num2str(i),' * Wd'])