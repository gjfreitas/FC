%% Problema 3.1: Oscilador harmónico — Runge–Kutta de 2ª ordem
%%
clear all
close all
clc

% a)
% w = (K/m)^1/2
% dv/dt = -K/m * x = -w^2 * x
% dx/dt = v

% b) Pela tabela:
%𝑟1 = 𝑓(t(k), y(k))
%𝑟2 = 𝑓 (t(k) + h/2, y(k) + 𝑟1 * h/2)
%y(k+1) = y(k) + r2 * h

% c)

h = 0.1;
x0 = 1;v0 = 0;K = 16;m = 1;
t0=0;tf=10;

w = sqrt(K/m);
t = t0:h:tf;
N = length(t);
x = nan(1,N);
v = nan(1,N);
x(1) = x0;v(1) = v0;

for k=1:N-1
    r1x = v(k);
    r1v = -w^2 * x(k);
    
    r2x = v(k) + r1v*h/2;
    r2v = -w^2 * ( x(k) + r1x*h/2 );
    
    v(k+1) = v(k) + r2v*h;
    x(k+1) = x(k) + r2x*h;
end

% Solução analitica 
xsa = x0*cos(w.*t); 
vsa = -w*x0*sin(w.*t);

figure(1)
subplot(1,2,1)
plot(t,x,'r--',t, xsa,'k.-'), xlabel('t'), ylabel('x')
subplot(1,2,2)
plot(t,v,'r--',t, vsa,'k.-'),xlabel('t'), ylabel('v')


%% d)

h = 0.1;
x0 = 1;v0 = 0;K = 16;m = 1;
t0=0;tf=10;

w = sqrt(K/m);
t = t0:h:tf;
N = length(t);
x = nan(1,N);
v = nan(1,N);
x(1) = x0;v(1) = v0;

fx = @(t,x,v) v; %derivada posicao
fv = @(t,x,v) -K*x/m; %derivada tempo
for k=1:N-1
    r1x = fx( t(k),x(k),v(k) );
    r1v = fv( t(k),x(k),v(k) );
    
    r2x = fx( t(k) +h/2 ,x(k) + r1x * h/2 ,v(k) + r1v*h/2 );
    r2v = fv( t(k) +h/2 ,x(k) + r1x * h/2 ,v(k) + r1v*h/2 );
    
    v(k+1) = v(k) + r2v*h;
    x(k+1) = x(k) + r2x*h;
end

% Solução analitica 
xsa = x0*cos(w.*t); 
vsa = -w*x0*sin(w.*t);

figure(1)
subplot(1,2,1)
plot(t,x,'r--',t, xsa,'k.-'), xlabel('t'), ylabel('x')
subplot(1,2,2)
plot(t,v,'r--',t, vsa,'k.-'),xlabel('t'), ylabel('v')

%% e)
clear all
close all
clc

h = 0.1;
x0 = 1;v0 = 0;K = 16;m = 1;
t0=0;tf=10;

w = sqrt(K/m);
t = t0:h:tf;
N = length(t);
x = nan(1,N);
v = nan(1,N);
x(1) = x0;v(1) = v0;

fx = @(V) V;
fv = @(X) -K*X/m;

for k=1:N-1
     
    r1x = fx(v(k));
    r1v = fv(x(k));
    r2x = fx(v(k)+r1v*h/2); 
    r2v = fv(x(k)+r1x*h/2);
    
    v(k+1) = v(k) + r2v*h;
    x(k+1) = x(k) + r2x*h;
end

% Solução analitica 
xsa = x0*cos(w.*t); 
vsa = -w*x0*sin(w.*t);

figure(1)
subplot(1,2,1)
plot(t,x,'r--',t, xsa,'k.-'), xlabel('t'), ylabel('x')
subplot(1,2,2)
plot(t,v,'r--',t, vsa,'k.-'),xlabel('t'), ylabel('v')


%% f)

clear all
close all
clc

h = 0.01;
x0 = 1;v0 = 0;K = 16;m = 1;
t0=0;tf=10;

w = sqrt(K/m);
t = t0:h:tf;
N = length(t);
x = nan(1,N);
v = nan(1,N);
x(1) = x0;v(1) = v0;

fx = @(v) v; %derivada posicao
fv = @(x) -K*x/m; %derivada tempo
for k=1:N-1
    r1x = fx( v(k) );
    r1v = fv( x(k) );
    
    r2x = fx( v(k) + r1v * h/2 );
    r2v = fv( x(k) + r1x * h/2 );
    
    v(k+1) = v(k) + r2v * h;
    x(k+1) = x(k) + r2x * h;
end

%euler
x_Euler = nan(1,N);
v_Euler = nan(1,N);
x_Euler(1) = x0;
v_Euler(1) = v0;
for k = 1:(N-1)%Euler
    v_Euler(k+1) = v_Euler(k)-(K/m)*x_Euler(k)*h;
    x_Euler(k+1) = x_Euler(k)+v_Euler(k)*h;
end

subplot(2,2,1)
plot( t,x,'.m',t, x0*cos(w.*t),'k.-',t,x_Euler,'g.'), xlabel('t'), ylabel('x')
subplot(2,2,2)
plot( t,v,'.m',t, -w*x0*sin(w.*t),'k.-',t,v_Euler,'g.'), xlabel('t'), ylabel('v')

subplot(2,2,3)
plot(x,v,x_Euler,v_Euler)

Ec = 1/2 * m *v.^2;
Ep = 1/2 * K *x.^2;
Em = Ec + Ep;

Ec_Euler = 1/2 * m *v_Euler.^2;
Ep_Euler = 1/2 * K *x_Euler.^2;
Em_Euler = Ec_Euler + Ep_Euler;

Ec_exact = 1/2 * m *( -w*x0*sin(w.*t) ).^2;
Ep_exact = 1/2 * K *(cos(w.*t)).^2;
Em_exact = Ec_exact + Ep_exact;

subplot(2,2,4)
plot(t,Em,t,Em_Euler,t,Em_exact)

