24%% Problema 3.2: Oscilador harmónico — Runge–Kutta de 4ª ordem
%% a)
clear all
close all
clc 

% Pela tabela
% r1 = f(x,v)
% r2 = f(x+r1*h/2,v+r1*h/2)
% r3 = f(x+r2*h/2,v+r2*h/2)
% r2 = f(x+r3*h,v+r3*h)
% v(k+1) = v(k) + 1/6 *(r1x + 2*r2x + 2*r3x +r4x)*h;
% x(k+1) = x(k) + 1/6 *(r1v + 2*r2v + 2*r3v +r4v)*h;

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

for k=1:(N-1)
     
    r1x = fx(v(k));
    r1v = fv(x(k));
    
    r2x = fx(v(k)+r1v*h/2); 
    r2v = fv(x(k)+r1x*h/2);
        
    r3x = fx(v(k)+r2v*h/2); 
    r3v = fv(x(k)+r2x*h/2);
    
    r4x = fx(v(k)+r3v*h); 
    r4v = fv(x(k)+r3x*h);
    
    x(k+1) = x(k) + 1/6 *(r1x + 2*r2x + 2*r3x +r4x)*h;
    v(k+1) = v(k) + 1/6 *(r1v + 2*r2v + 2*r3v +r4v)*h;
end
figure(1)
subplot(1,2,1)
plot( t,x,t, x0*cos(w.*t),'k.-'), xlabel('t'), ylabel('x')
subplot(1,2,2)
plot( t,v,t, -w*x0*sin(w.*t),'k.-'), xlabel('t'), ylabel('v')



Ec = 1/2 * m *v.^2;
Ep = 1/2 * K *x.^2;
Em_RK4 = Ec + Ep;

Ec_exact = 1/2 * m *( -w*x0*sin(w.*t) ).^2;
Ep_exact = 1/2 * K *(cos(w.*t)).^2;
Em_exact = Ec_exact + Ep_exact;

figure(2)
plot(t,Em_RK4,t,Em_exact)

%% d)
clear all
close all
clc 

% Pela tabela
% r1 = f(x,v)
% r2 = f(x+r1*h/2,v+r1*h/2)
% r3 = f(x+r2*h/2,v+r2*h/2)
% r2 = f(x+r3*h,v+r3*h)
% v(k+1) = v(k) + 1/6 *(r1x + 2*r2x + 2*r3x +r4x)*h;
% x(k+1) = x(k) + 1/6 *(r1v + 2*r2v + 2*r3v +r4v)*h;


x0 = 1;v0 = 0;K = 16;m = 1;
t0=0;tf=10;
w = sqrt(K/m);
h = [2E-1, 1E-1,5E-2, 2E-2,1E-2, 5E-3,2E-3,1E-3]; %s
% h = [1E-2, 5E-3,1E-4, 5E-4,1E-4, 5E-5,1E-5]; %s
erro = nan(length(h),1);
x_exact_t10 = x0*cos(w*10); % solução exata para t = 10s
% v_exact_t10 = -w*x0*sin(w*10);

for index_h = 1:length(h)
    t = t0:h(index_h):tf;
    
    N = length(t);
    x = nan(N,1);
    v = nan(N,1);
    x(1) = x0;v(1) = v0;

    fx = @(V) V;
    fv = @(X) -K*X/m;

    for k=1:(N-1)

        r1x = fx(v(k));
        r1v = fv(x(k));

        r2x = fx(v(k)+r1v*h(index_h)/2); 
        r2v = fv(x(k)+r1x*h(index_h)/2);

        r3x = fx(v(k)+r2v*h(index_h)/2); 
        r3v = fv(x(k)+r2x*h(index_h)/2);

        r4x = fx(v(k)+r3v*h(index_h)); 
        r4v = fv(x(k)+r3x*h(index_h));

        x(k+1) = x(k) + 1/6 *(r1x + 2*r2x + 2*r3x +r4x)*h(index_h);
        v(k+1) = v(k) + 1/6 *(r1v + 2*r2v + 2*r3v +r4v)*h(index_h);
    end
    erro(index_h) = abs(x_exact_t10 - x(N));
    
end


plot(log10(h),log10(erro),'ko'), lsline
lsline %traçar uma recta sobre a curva que obteve
aux = polyfit(log10(h),log10(erro),1)  %1 corresponde a ordem do polinomio, como é uma reta a ordem = 1
declive = aux(1);
disp([' declive: ',num2str(declive)])
