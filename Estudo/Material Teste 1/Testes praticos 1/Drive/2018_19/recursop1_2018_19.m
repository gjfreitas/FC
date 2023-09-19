%% Exame Prático de Recurso — Parte 1 Física Computacional — 2018/2019 9 de julho de 2019

% 1 
clear all
close all 
clc

R = 0.2; %m
g = 9.8; %m s^2

% 1 a)
w = 5; %rad/s
theta0 = pi/4; %rad
dtheta0 = 0;

t0 = 0;
tf = 10;
h = 0.01;
t = t0:h:tf;
N = length(t);

theta = zeros(1,N); theta(1) = theta0;
dtheta = zeros(1,N); dtheta(1) = dtheta0;

ftheta = @(dtheta) dtheta;            
fdtheta = @(theta) (w^2*cos(theta)-g/R)*sin(theta);

% Runge Kutta 4

for i = 1:N-1
    r1dt = fdtheta(theta(i));                 
    r1t = ftheta(dtheta(i));  
    
    r2dt = fdtheta(theta(i) + r1t * h/2);     
    r2t = ftheta(dtheta(i) + r1dt * h/2);
    
    r3dt = fdtheta(theta(i)+ r2t * h/2);
    r3t = ftheta(dtheta(i) + r2dt * h/2);
    
    r4dt = fdtheta(theta(i) + r3t * h);
    r4t = ftheta(dtheta(i) + r3dt * h);
    
    theta(i+1) = theta(i) + 1/6*(r1t + 2*r2t + 2*r3t + r4t)*h;
    dtheta(i+1) = dtheta(i) + 1/6*(r1dt + 2*r2dt + 2*r3dt + r4dt)*h;
end

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

plot(t, theta), xlabel('tempo (s)'), ylabel('theta (rad)')
title('Método Runge-Kutta 4a ordem')

%% b)
clear all
close all 
clc

R = 0.2; %m
g = 9.8; %m s^2

% 1 a)
w = [5, 8, 9]; %rad/s
theta0 = pi/4; %rad
dtheta0 = 0;

for k = 1: length(w)
    t0 = 0;
    tf = 10;
    h = 0.01;
    t = t0:h:tf;
    N = length(t);

    theta = zeros(1,N); theta(1) = theta0;
    dtheta = zeros(1,N); dtheta(1) = dtheta0;

    ftheta = @(dtheta) dtheta;            
    fdtheta = @(theta) (w(k)^2*cos(theta)-g/R)*sin(theta);

% Runge Kutta 4

    for i = 1:N-1
        r1dt = fdtheta(theta(i));                 
        r1t = ftheta(dtheta(i));  

        r2dt = fdtheta(theta(i) + r1t * h/2);     
        r2t = ftheta(dtheta(i) + r1dt * h/2);

        r3dt = fdtheta(theta(i)+ r2t * h/2);
        r3t = ftheta(dtheta(i) + r2dt * h/2);

        r4dt = fdtheta(theta(i) + r3t * h);
        r4t = ftheta(dtheta(i) + r3dt * h);

        theta(i+1) = theta(i) + 1/6*(r1t + 2*r2t + 2*r3t + r4t)*h;
        dtheta(i+1) = dtheta(i) + 1/6*(r1dt + 2*r2dt + 2*r3dt + r4dt)*h;
    end

        figure(1)
        hold on

        if w(k) == 5
            plot(t,theta,'b-'), xlabel('t(s)'), ylabel('theta')
        end

        if w(k) == 8
            plot(t,theta,'m-'), xlabel('t(s)'), ylabel('theta')
        end

        if w(k) == 9
            plot(t,theta,'k-'), xlabel('t(s)'), ylabel('theta')
        end

end
legend('w = 5','w = 8','w = 9'), xlabel('tempo (s)'), ylabel('theta (rad)')
title('Método Runge-Kutta 4a ordem')

%% 1 c)
clear all
close all 
clc

R = 0.2; %m
g = 9.8; %m s^2

% 1 a)
w = 9; %rad/s
theta0 = acos(g/(w^2 * R)); %rad
dtheta0 = 0;

t0 = 0;
tf = 10;
h = 0.01;
t = t0:h:tf;
N = length(t);

theta = zeros(1,N); theta(1) = theta0;
dtheta = zeros(1,N); dtheta(1) = dtheta0;

ftheta = @(dtheta) dtheta;            
fdtheta = @(theta) (w^2*cos(theta)-g/R)*sin(theta);

% Runge Kutta 4

for i = 1:N-1
    r1dt = fdtheta(theta(i));                 
    r1t = ftheta(dtheta(i));  
    
    r2dt = fdtheta(theta(i) + r1t * h/2);     
    r2t = ftheta(dtheta(i) + r1dt * h/2);
    
    r3dt = fdtheta(theta(i)+ r2t * h/2);
    r3t = ftheta(dtheta(i) + r2dt * h/2);
    
    r4dt = fdtheta(theta(i) + r3t * h);
    r4t = ftheta(dtheta(i) + r3dt * h);
    
    theta(i+1) = theta(i) + 1/6*(r1t + 2*r2t + 2*r3t + r4t)*h;
    dtheta(i+1) = dtheta(i) + 1/6*(r1dt + 2*r2dt + 2*r3dt + r4dt)*h;
end

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

plot(t, theta), xlabel('tempo (s)'), ylabel('theta (rad)')
title('Método Runge-Kutta 4a ordem')

%% 1 d)
clear all
close all 
clc

R = 0.2; %m
g = 9.8; %m s^2

% 1 a)
w = 9; %rad/s
theta0 = acos(g/(w^2 * R)); %rad
dtheta0 = 0;

t0 = 0;
tf = 10;
h = 0.01;
t = t0:h:tf;
N = length(t);


% ode45
reltol = 3E-14;
abstol_1 = 1E-13;
abstol_2 = 1E-13;
options = odeset('RelTol',reltol,'AbsTol',[abstol_1 abstol_2]);
[t_ode45,sol] = ode45(@f,[t0 tf],[theta0 dtheta0],options);

figure(1)
plot(t_ode45,sol(:,1),'r--'), xlabel('t'), ylabel('x')
figure(2)
plot(t_ode45,sol(:,2),'r--'), xlabel('t'), ylabel('v')

%% 2 a)
clear all
close all
clc

lbd = 2.2;

xi = 0;
xf = pi;
h = 0.01;
x = xi:h:xf;
n = length(x);


y = zeros(1,n); y(1) = 0;
vy = zeros(1,n); vy(1) = 1;

% Euler-Cromer
%dy/dx = vy -> fy = vy

fy = @(vy) vy;
fvy = @ (x, y, vy) ((x-lbd)*y+1.6*sin(x)*cos(x)*vy)/(1-0.8*sin(x)^2);


for i = 1:n-1
    vy(i+1) = vy(i) + fvy( x(i), y(i), vy(i)) * h;
    y(i+1) = y(i) + vy(i+1) * h;            %metodo de euler-cromer usa vx(i+1)
end



%grafico y(x)
plot(x,y, 'm'), xlabel('x(m)'), ylabel('y(m)')


%% 2 b) e 2 c) Não sai para o teste 1 2020/21
