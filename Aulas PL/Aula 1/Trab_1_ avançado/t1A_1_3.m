%% Trabalho 1- Avançado(2020/21)
%% 1.3

clear all
close all
clc

% Constantes
g = 9.81;
m = 450*10^(-3); %kg 
r = (70/(2*pi))*10^(-2); %m
theta = (10*pi)/180; 
rho = 1.225; %kg/m^3
A = pi*r^2;
w = (600*2*pi)/60;
Cm = 1;
B = (1/2)*Cm*rho*A*r;


h = 0.001;
t = 0:h:10;

N = length(t);
vx = zeros(1,N);
vy = zeros(1,N);
vz = zeros(1,N);
v = zeros(1,N);
v(1) = (80 * 10^3)/3600; %m/s
vx(1) = v(1)*cos(theta);
vz(1) = v(1)*sin(theta);


x = zeros(1,N);
y = zeros(1,N);
z = zeros(1,N);

% Em cada direção, atuam as seguintes forças:
%Fx = Fdx + FLx; Fy = Fly + Fdy; Fz = Fg + Fdz

% Método de Euler

for i = 1:(length(t)-1)
    if abs(v(i)) <= 9
        ax = (-B*w*vy(i)-(0.015*abs(v(i))*vx(i)))/m;
        vx(i+1) = vx(i) + ax*h;
       
        ay =  (-(0.015*abs(v(i))*vy(i)) + B*w*vx(i))/m;
        vy(i+1) = vy(i) + ay*h;
        
        az = (-(0.015*abs(v(i))*vz(i)))/m - g;
        vz(i+1) = vz(i) + az*h;
        
    elseif abs(v(i)) <= 20
        ax = ((-(0.25147 + 0.17431*abs(v(i)) - 0.01384*abs(v(i))^2 + 0.00054*abs(v(i))^3)*vx(i))/abs(v(i)) - B*w*vy(i))/m;
        vx(i+1) = vx(i) + ax*h;
        
        ay = ((-(0.25147 + 0.17431*abs(v(i)) - 0.01384*abs(v(i))^2 + 0.00054*abs(v(i))^3)*vy(i))/abs(v(i)) + B*w*vx(i))/m;
        vy(i+1) = vy(i) + ay*h;
        
        az = (-(0.25147 + 0.17431*abs(v(i)) - 0.01384*abs(v(i))^2 + 0.00054*abs(v(i))^3)*vz(i))/(m*abs(v(i))) - g;
        vz(i+1) = vz(i) + az*h;
    
    else
        ax = ((-(-4.025 + 0.323*abs(v(i)))*vx(i))/abs(v(i)) - B*w*vy(i))/m;
        vx(i+1) = vx(i) + ax*h;
        
        ay = ((-(-4.025 + 0.323*abs(v(i)))*vy(i))/abs(v(i)) + B*w*vx(i))/m;
        vy(i+1) = vy(i) + ay*h;
        
        az =(-(-4.025 + 0.323*abs(v(i)))*vz(i))/(m*abs(v(i))) - g;
        vz(i+1) = vz(i) + az*h;
        
    end
    
    x(i+1) = x(i) + vx(i)*h;
    y(i+1) = y(i) + vy(i)*h;
    z(i+1) = z(i) + vz(i)*h;
    
    v(i+1)=norm(vx(i+1)+vy(i+1)+ vz(i+1));
   
    if z(i+1) < 0
        break
    end
end

% plot3(x,y,z,'o')

%inicialmente assumiu-se vx(1) = v*cos(THETA) e vy(1) = 0. Portanto, se não
%houver desvio lateral espera-se quer vy(tf) = 0

vx_f = vx(i);
vy_f = vy(i)
PHI = atan(vy_f/vx_f);

xmax = x(i) % distancia percorrida
dlat = y(i) % desvio lateral