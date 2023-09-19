%% %% Exame Prático — Física Computacional — 2018/2019 21 de junho de 2019

% 1 a)
clear all
close all
clc

vlim = 6.8;
g = 9.8;
ti = 0;
tf = 10;
v0 = 30;
v0x = v0 * cosd(60);
v0z = v0 * sind(60);
x0 = 0;
z0 = 1;

h = 0.01;
t = ti:h:tf;
N = length(t);
vx = zeros(1,N);
vz = zeros(1,N);
z = zeros(1,N);
x = zeros(1,N);

x(1) = x0;
z(1) = z0;
vx(1) = v0x;
vz(1) = v0z;

fz = @(vz) vz;            %derivada da altura -> velocidade em z
fvz = @(vz) (- g/(vlim)^2 * vz*abs(vz) - g);       %derivada da velocidade em z
fx = @(vx) vx;
fvx = @(vx) (- g/(vlim)^2 * vx*abs(vx));

% Runge Kutta 2 ordem

for k = 1:N-1
    r1vx = fvx(vx(k));  
    r1vz = fvz(vz(k));
    r1z = fz(vz(k));
    r1x = fx(vx(k));
    
    r2vx = fvx(vx(k)+r1vx*h);  
    r2vz = fvz(vz(k)+r1vz*h);
    r2z = fz(vz(k)+r1vz*h);
    r2x = fx(vx(k)+r1vx*h);
    
    vx(k+1) = vx(k) + r1vx*h/2 + r2vx*h/2;
    vz(k+1) = vz(k) + r2vx*h/2 + r2vz*h/2;
    x(k+1) = x(k) + r1x*h/2 + r2x*h/2;
    z(k+1) = z(k) + r1z*h/2 + r2z*h/2;

    if z(k) < 0
        break
    end
end 

%interpolacao para descobrir o ponto em que a função z(t) é 0
t_impacto = interp1([z(k-1),z(k)],[t(k-1),t(k)],0);
x_impacto = interp1([t(k-1),t(k)],[x(k-1),x(k)],t_impacto);
disp(['tempo de voo: ',num2str(t_impacto),' s'])
disp(['Alcance : ',num2str(x_impacto),' m'])

figure(1)
plot(x(1:k),z(1:k))
title('Trajetória do volante')
xlabel('x (m)')
ylabel('altura (m)')

figure(2)
plot(t(1:k),z(1:k),'r')
title('Evolução da altura do volante')
xlabel('tempo (s)')
ylabel('altura (m)')

%% 1 b)
clear all
close all
clc

vlim = 6.8;
g = 9.8;
ti = 0;
tf = 2.5;
v0 = 30;
v0x = v0 * cosd(60);
v0z = v0 * sind(60);
x0 = 0;
z0 = 1;
h = 0.01;
t = ti:h:tf;


reltol = 3E-14;
abstol_1=1E-13;
abstol_2=1E-13;


options = odeset('RelTol', reltol, 'AbsTol', [abstol_1 abstol_2 abstol_1 abstol_2]);      %opcoes para o ode45
[t, sol] = ode45(@f, t, [x0 v0x z0 v0z], options); 

plot(sol(:,1), sol(:,3), '-k')
title('Volante de badminton - ode45')
xlabel('Abcissa (m)')
ylabel('Altura (m)')



