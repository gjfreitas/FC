%% Exame Prático Física Computacional — 2015/2016
%% 2 a)
clear all
close all
clc

L = 50;
xf = L;
x0 = 0;

c = 0.094;
k = 0.93;
p = 8.9;
a = 0.02;
Tamb = 20;
const = (a*c*p)/k;

T0 = 10;
% dT/dx = Q
Q0 = 2;

%ODE45
reltol = 1*10^(-10);
abstol_1 = 1*10^(-10);
abstol_2 = 1*10^(-10);
options = odeset('RelTol',reltol,'AbsTol',[abstol_1 abstol_2]);


[x,sol] = ode45(@f,[x0 xf],[T0 Q0],options,const);
T = sol(:,1);
plot(x,T), xlabel('x(m)'), ylabel('T(ºC ) ')

%% 2 b)
 clear all
 close all
 clc
 
L = 50;
c = 0.094;
k = 0.93;
p = 8.9;
a = 0.02;
Tamb = 20;
const = (a*c*p)/k;
 
 x0 = 0;
 xf = L;
 h = 0.01;
 x = x0:h:xf;
 N = length(x);
 T0 = 200;
%dT/dx = Q
 
 T = nan(1,N); T(1) = T0;
 Q = nan(1,N);
 
 % Shooting
 guess = [0 1.2];
 B = 50;
 nshots = 10000;
 
 for ishot = 1:nshots
     Q = nan(1,N);
     Q(1) = guess(ishot);
     for k = 1:N-1
         Q(k+1) = Q(k) + (const*(T(k)-Tamb))*h;
         T(k+1) = T(k) + Q(k+1) *h;
     end
     result(ishot) = T(end);
     
     if ishot >= 2
        m = (result(ishot)-result(ishot-1))/(guess(ishot)-guess(ishot-1));
        guess(ishot+1) = guess(ishot) + (B - result(ishot))/m;
        if abs(guess(ishot) - guess(ishot-1)) < 1E-6
            break
        end   
     end
 end
 
solucao = guess(ishot);
disp(['dT/dx(0) => T(L) = ',num2str(result(ishot)),' : ',num2str(guess(ishot))])
 plot(x,T)
 
 % 2 c)
 
%ODE45
reltol = 1*10^(-10);
abstol_1 = 1*10^(-10);
abstol_2 = 1*10^(-10);
options = odeset('RelTol',reltol,'AbsTol',[abstol_1 abstol_2]);
 
 Q0 = solucao;
 [x,sol] = ode45(@f,[x0 xf],[T0 Q0],options,const);
 dx = (xf-x0)/(length(x)-1);
 T = sol(:,1);
 Q = sol(:,2);
 
 troca_int = k .* Q;
 troca_ext = c*p*a*dx.*(T - Tamb);
 
 plot(x,troca_int,'r.-',x,troca_ext,'b.-') % não sao simetricas ...