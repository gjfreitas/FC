%% Trabalho 1- Avançado(2020/21)
%% 1.2

clear all
close all
clc

g = 9.81;
m = 57*10^(-3); %kg 
R = (67/2)*10^(-3); %m
theta = (5*pi)/180; 
rho = 1.225; %kg/m^3
A = pi*R^2;
B = (1/2)*rho*A;

t = 0:0.001:10;
N = length(t);
vx = zeros(1,N);
vz = zeros(1,N);
v = zeros(1,N);

x = zeros(1,N);
z = zeros(1,N);
v(1) = 20; %m/s
vx(1) = v(1)*cos(theta);
vz(1) = v(1)*sin(theta);
z(1) = 0.7; %m


% % alínea a - sem rotação
% w = 0;
% signal_w = 1;
% % alínea b - topspin
w = 3000*2*pi/60;
signal_w = 1;
% % alínea c - backspin
% w = -3000*2*pi/60;
% signal_w = -1;

% Método de Euler
h = 0.001;
for i = 1:(N-1)
        S = (R*abs(w))/abs(v(i));
        Cd = 0.508 + (22.503 + 4.196*(S^(-2.5)))^(-0.4);
        Cl = (2.022 + 0.981*S^(-1))^(-1);
        
        ax = (B/m)*abs(v(i))*(Cl*vz(i)*signal_w-Cd*vx(i));
        vx(i+1) = vx(i) + ax*h;
        x(i+1) = x(i) + vx(i)*h;
        
        az = (-B/m)*abs(v(i))*(Cl*vx(i)*signal_w+Cd*vz(i))- g;
        vz(i+1)= vz(i) + az*h;
        z(i+1) = z(i) + vz(i)*h;
        
        v(i+1)=sqrt(vx(i+1)^2 + vx(i+1)^2);
        
        if z(i+1) < 0
            break
        end
end

plot(x,z,'o')

% interpolação de Lagrange
[z_max,ind]=max(z);
aux=lagr(x(ind-1:ind+1),z(ind-1:ind+1));

