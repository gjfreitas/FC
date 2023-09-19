%% Problema 4.2: Oscilador de van der Pol — Runge–Kutta de 4ª ordem
%% a1) F0 = 1.0
clear all
close all
clc

t0 = 0; %s
h = 0.01; %s
tf = 100; %s
v0 = 2; %m/s
y0 = 0; %m
eps = 0.1;
F0 = 1.0;

t = t0:h:tf;
N = length(t);
y = nan(1,N);
v = nan(1,N);
y(1) = y0;v(1) = v0;

fy = @(V) V;
fv = @(T,Y,V) F0*cos(1.7*T)-eps.*(Y.^2-1).*V-Y;

for k=1:(N-1)
    r1y = fy(v(k));
    r1v = fv(t(k),y(k),v(k));
    
    r2y = fy(v(k) + r1v * h/3);
    r2v = fv(t(k)+h/3,y(k) + r1y * h/3,v(k) + r1v * h/3);
        
    r3y = fy(v(k) - r1v * h/3 + r2v * h);
    r3v = fv(t(k) + h* 2/3,y(k) - r1y * h/3 + r2y * h,v(k) - r1v * h/3 + r2v * h);
    
    r4y = fy(v(k) + r1v * h - r2v * h + r3v *h);
    r4v = fv(t(k) + h,y(k) + r1y * h - r2y * h + r3y * h,v(k) + r1v * h - r2v * h + r3v *h);
    
    y(k+1) = y(k) + 1/8 *(r1y + 3 * r2y + 3 * r3y + r4y) * h;
    v(k+1) = v(k) + 1/8 *(r1v + 3 * r2v + 3 * r3v + r4v) * h;
end

% figure(1)
% plot(t,y,'r-'), xlabel('t'), ylabel('y')
% figure(2)
% plot(y,v,'r-'),xlabel('y'), ylabel('v')
% figure(3)
% plot(t,v,'r-'),xlabel('t'), ylabel('v')

figure(1)
subplot(1,3,1)
plot(t,y,'r-'), xlabel('t'), ylabel('y')
subplot(1,3,2)
plot(t,v,'r-'),xlabel('t'), ylabel('v')
subplot(1,3,3)
plot(y,v,'r-'),xlabel('y'), ylabel('v')

%% a1) F0 = 1.5
clear all
close all
clc

t0 = 0; %s
h = 0.01; %s
tf = 100; %s
v0 = 2; %m/s
y0 = 0; %m
eps = 0.1;
F0 = 1.5;

t = t0:h:tf;
N = length(t);
y = nan(1,N);
v = nan(1,N);
y(1) = y0;v(1) = v0;

fy = @(V) V;
fv = @(T,Y,V) F0*cos(1.7*T)-eps.*(Y.^2-1).*V-Y;

for k=1:(N-1)
    r1y = fy(v(k));
    r1v = fv(t(k), y(k), v(k));
    
    r2y = fy(v(k) + r1v * h/3);
    r2v = fv(t(k)+ h/3, y(k) + r1y * h/3, v(k) + r1v * h/3);
        
    r3y = fy(v(k) - r1v * h/3 + r2v * h);
    r3v = fv(t(k) + h * 2/3, y(k) - r1y * h/3 + r2y * h, v(k) - r1v * h/3 + r2v * h);
    
    r4y = fy(v(k) + r1v * h - r2v * h + r3v *h);
    r4v = fv(t(k) + h, y(k) + r1y * h - r2y * h + r3y * h, v(k) + r1v * h - r2v * h + r3v *h);
    
    y(k+1) = y(k) + 1/8 *(r1y + 3 * r2y + 3 * r3y + r4y) * h;
    v(k+1) = v(k) + 1/8 *(r1v + 3 * r2v + 3 * r3v + r4v) * h;
end

% figure(1)
% plot(t,y,'r-'), xlabel('t'), ylabel('y')
% figure(2)
% plot(y,v,'r-'),xlabel('y'), ylabel('v')
% figure(3)
% plot(t,v,'r-'),xlabel('t'), ylabel('v')

figure(1)
subplot(1,3,1)
plot(t,y,'r-'), xlabel('t'), ylabel('y')
subplot(1,3,2)
plot(t,v,'r-'),xlabel('t'), ylabel('v')
subplot(1,3,3)
plot(y,v,'r-'),xlabel('y'), ylabel('v')