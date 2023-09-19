%% Exame Prático de Recurso Física Computacional — 2015/2016
%% 2 a)
clear all
close all
clc

% Constantes
m = 145E-3; %kg
R = 36.3E-3; %m
y0 = 0.9;
x0 = 0;
ang = 30;
v0 = 45;
po = 1.225;
delta = 5;
vdelta = 35;
g = 9.8;

v0x = v0*cosd(ang);
v0y = v0*sind(ang);
v = [v0x v0y];

 t0 = 0; %s
 tf = 10; %s
 h = 0.01; %s
 t = t0:h:tf; %s
 
 N = length(t);
 x = nan(1,N); x(1) = x0;
 y = nan(1,N); y(1) = y0;
 vx = nan(1,N); vx(1) = v0x;
 vy = nan(1,N);vy(1) = v0y;

 A = pi*R^2;
 const = 1/2 * po * A;
 
 
% g so atua em relação a y
for k = 1:N-1
     if y(k) < 0
         break
     end
    
    Cd(k) = 0.2194 + 0.3263/(1 + exp((sqrt(vy(k)^2 + vx(k)^2) - vdelta)/delta));
    
    vx(k+1) = vx(k) + (-const * Cd(k) * vx(k))*h;
    x(k+1) = x(k) + vx(k+1)*h;

    vy(k+1) = vy(k) + (-const * Cd(k) * vy(k) - g)*h;
    y(k+1) = y(k)+ vy(k+1)*h;

end
 
plot(x,y,'m-'),xlabel('t'), ylabel('y'), title('Trajetória')

%% 2 b)
clear all
close all
clc

% Constantes
m = 145E-3; %kg
R = 36.3E-3; %m
y0 = 0.9;
x0 = 0;
ang = 30;
po = 1.225;
delta = 5;
vdelta = 35;
g = 9.8;

A = pi*R^2;
const = 1/2 * po * A;

t0 = 0; %s
tf = 10; %s
h = 0.01; %s
t = t0:h:tf; %s

N = length(t);
x = nan(1,N); x(1) = x0;
y = nan(1,N); y(1) = y0;

 
% shooting
B = 100;
guess = [30.0 30.5];
tol = 1E-5;
nshots = 1000;

for j = 1:nshots
    clear v0
    clear v0y
    clear v0x
    
    v0 = guess(j);
    v0x = v0*cosd(ang);
    v0y = v0*sind(ang);
    vx = nan(1,N); vx(1) = v0x;
    vy = nan(1,N);vy(1) = v0y;

    for k = 1:N-1
        if y(k) < 0
            break
        end
        vx(k+1) = vx(k) + (-const*sqrt((vy(k)^2 + vx(k)^2)) * vx(k))*h;
        x(k+1) = x(k) + vx(k+1)*h;

        vy(k+1) = vy(k) + (-const * sqrt((vy(k)^2 + vx(k)^2)) * vy(k) - g)*h;
        y(k+1) = y(k)+ vy(k+1)*h;
    end

    result(j) = x(k);
    diff = B - result(j);

    if j>= 2
        m = (result(j)-result(j-1))/(guess(j)-guess(j-1));
        guess(j+1) = guess(j)+(diff)/m;
        if abs(guess(j)-guess(j-1)) < tol
                break
        end
    end
  
end

disp(['alcance: ',num2str(result(j)),' m'])
disp(['|v0|: ',num2str(guess(j)),' m/s'])
plot(x,y,'m-'),xlabel('x'), ylabel('y'), title('Trajetória')