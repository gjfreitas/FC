%% Problema 4.1: Oscilador de van der Pol — Runge–Kutta de 3ª ordem

%% a1) e = 0.1
clear all
close all
clc

t0 = 0; %s
h = 0.01; %s
tf = 100; %s
v0 = 0.7; %m/s
y0 = 0.2; %m
eps = 0.1;

t = t0:h:tf;
N = length(t);
y = nan(1,N);
v = nan(1,N);
y(1) = y0;v(1) = v0;

fy = @(V) V;
fv = @(Y,V) -eps.*(Y.^2-1).*V-Y;

for k=1:(N-1)
    r1y = fy(v(k));
    r1v = fv(y(k),v(k));
    
    r2y = fy(v(k) + r1v * h);
    r2v = fv(y(k) + r1y * h,v(k) + r1v * h);
        
    r3y = fy(v(k) + r1v * h/4 + r2v * h/4);
    r3v = fv(y(k) + r1y * h/4 + r2y * h/4,v(k) + r1v * h/4 + r2v * h/4);
    
    y(k+1) = y(k) + 1/6 *(r1y + r2y + 4 * r3y) * h;
    v(k+1) = v(k) + 1/6 *(r1v + r2v + 4 * r3v) * h;
end

figure(1)
plot(t,y,'b-'), xlabel('t'), ylabel('y')
figure(2)
plot(y,v,'b-'),xlabel('y'), ylabel('v')

%% a2) e = 1
clear all
close all
clc

t0 = 0; %s
h = 0.01; %s
tf = 100; %s
v0 = 0.7; %m/s
y0 = 0.2; %m
eps = 1;

t = t0:h:tf;
N = length(t);
y = nan(1,N);
v = nan(1,N);
y(1) = y0;v(1) = v0;

fy = @(V) V;
fv = @(Y,V) -eps.*(Y.^2-1).*V-Y;

for k=1:(N-1)
    r1y = fy(v(k));
    r1v = fv(y(k),v(k));
    
    r2y = fy(v(k) + r1v * h);
    r2v = fv(y(k) + r1y * h,v(k) + r1v * h);
        
    r3y = fy(v(k) + r1v * h/4 + r2v * h/4);
    r3v = fv(y(k) + r1y * h/4 + r2y * h/4,v(k) + r1v * h/4 + r2v * h/4);
    
    y(k+1) = y(k) + 1/6 *(r1y + r2y + 4 * r3y) * h;
    v(k+1) = v(k) + 1/6 *(r1v + r2v + 4 * r3v) * h;
end

figure(1)
plot(t,y,'b-'), xlabel('t'), ylabel('y')
figure(2)
plot(y,v,'b-'),xlabel('y'), ylabel('v')

%% b1) e = 0.1,v(0) = 7,y(o) = 2
clear all
close all
clc

t0 = 0; %s
h = 0.01; %s
tf = 100; %s
v0 = 7; %m/s
y0 = 2; %m
eps = 0.1;

t = t0:h:tf;
N = length(t);
y = nan(1,N);
v = nan(1,N);
y(1) = y0;v(1) = v0;

fy = @(V) V;
fv = @(Y,V) -eps.*(Y.^2-1).*V-Y;

for k=1:(N-1)
    r1y = fy(v(k));
    r1v = fv(y(k),v(k));
    
    r2y = fy(v(k) + r1v * h);
    r2v = fv(y(k) + r1y * h,v(k) + r1v * h);
        
    r3y = fy(v(k) + r1v * h/4 + r2v * h/4);
    r3v = fv(y(k) + r1y * h/4 + r2y * h/4,v(k) + r1v * h/4 + r2v * h/4);
    
    y(k+1) = y(k) + 1/6 *(r1y + r2y + 4 * r3y) * h;
    v(k+1) = v(k) + 1/6 *(r1v + r2v + 4 * r3v) * h;
end

figure(1)
plot(t,y,'m-'), xlabel('t'), ylabel('y')
figure(2)
plot(y,v,'m-'),xlabel('y'), ylabel('v')

%% b2) e = 1,v(0) = 7,y(o) = 2
clear all
close all
clc

t0 = 0; %s
h = 0.01; %s
tf = 100; %s
v0 = 7; %m/s
y0 = 2; %m
eps = 1;

t = t0:h:tf;
N = length(t);
y = nan(1,N);
v = nan(1,N);
y(1) = y0;v(1) = v0;

fy = @(V) V;
fv = @(Y,V) -eps.*(Y.^2-1).*V-Y;

for k=1:(N-1)
    r1y = fy(v(k));
    r1v = fv(y(k),v(k));
    
    r2y = fy(v(k) + r1v * h);
    r2v = fv(y(k) + r1y * h,v(k) + r1v * h);
        
    r3y = fy(v(k) + r1v * h/4 + r2v * h/4);
    r3v = fv(y(k) + r1y * h/4 + r2y * h/4,v(k) + r1v * h/4 + r2v * h/4);
    
    y(k+1) = y(k) + 1/6 *(r1y + r2y + 4 * r3y) * h;
    v(k+1) = v(k) + 1/6 *(r1v + r2v + 4 * r3v) * h;
end

figure(1)
plot(t,y,'m-'), xlabel('t'), ylabel('y')
figure(2)
plot(y,v,'m-'),xlabel('y'), ylabel('v')