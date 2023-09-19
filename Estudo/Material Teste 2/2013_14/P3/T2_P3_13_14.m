%% 2o Teste Prático Física Computacional — 2013/2014 —  Turma P3
%% 1 a)
 clear all
 close all
 clc
 
%  const = 0.1;
%  Tamb = 0;
 L = 5;
 
 x0 = 0;
 xf = L;
 T0 = 200;
 dT0 = -11.0;
 
%ODE45
reltol = 1*10^(-10);
abstol_1 = 1*10^(-10);
abstol_2 = 1*10^(-10);
options = odeset('RelTol',reltol,'AbsTol',[abstol_1 abstol_2]);


[x,sol] = ode45(@f,[x0 xf],[T0 dT0],options);
T = sol(:,1);
plot(x,T), xlabel('x(m)'), ylabel('T(ºC ) ')

%% 1 b)
 clear all
 close all
 clc
 
 const = 0.1;
 Tamb = 0;
 L = 5;
 
 x0 = 0;
 xf = L;
 h = 0.01;
 x = x0:h:xf;
 N = length(x);
 T0 = 200;
% dT/dx = Q
 
 T = nan(1,N); T(1) = T0;
 Q = nan(1,N);
 
 % Shooting
 guess = [-11.0 -10.2];
 B = 100;
 nshots = 100;
 
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
 
disp(['T(5) : ',num2str(result(ishot))])
disp(['Derivada de T(0) : ',num2str(guess(ishot))])
 plot(x,T)
 
 %% 2 a)
 clear all
 close all
 clc
 
 dt = 0.01;
 t = -10:dt:10;
 N = length(t);
 T = 0.26;
 
 B = zeros(1,N);
 
 for k = 1:N
     if t(k) > -T/2 && t(k) < T/2
         B(k) = 1;
     end
 end
 
 Bf = dt * fftshift(fft(B));
 
% 2 b) 
wmax = pi/dt;
dw = 2*wmax/N;
w = -N/2*dw:dw:(N/2-1)*dw;

Bf_exact = 2./w .* sin(T.*w/2);


plot(t,abs(Bf),'r.-',t,Bf_exact,'b.-')
 

 
 
 