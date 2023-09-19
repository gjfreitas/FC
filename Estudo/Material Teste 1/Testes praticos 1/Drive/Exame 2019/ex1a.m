% Física Computacional 2018/2019   Tiago Fernandes
% Problema FR 1.3 a)

clear all
clc

vlim = 6.8;
g = 9.8;
ti = 0;
tf = 10;

h = 0.01;
t = ti:h:tf;
n = length(t);
vx = zeros(1,n);
vz = zeros(1,n);
z = zeros(1,n);
x = zeros(1,n);
v0 = 30;
x(1) = 0;
z(1) = 1;
vx(1) = v0 * cosd(60);
vz(1) = v0 * sind(60);

fz = @(vz) vz;            %derivada da altura -> velocidade em z
fvz = @(vz) (- g/(vlim)^2 * vz*abs(vz) - g);       %derivada da velocidade em z
fx = @(vx) vx;
fvx = @(vx) (- g/(vlim)^2 * vx*abs(vx));

% metodo de runge kutta 2a ordem
for i = 1:n-1
    r1vx = fvx(vx(i));                 %sempre igual no RK
    r1vz = fvz(vz(i));
    r1z = fz(vz(i));                 %sempre igual no RK
    r1x = fx(vx(i));
    r2vx = fvx(vx(i)+r1vx*h);         %x somam aos rx
    r2vz = fvz(vz(i)+r1vz*h);
    r2z = fz(vz(i)+r1vz*h);         % v somam aos rv
    r2x = fx(vx(i)+r1vx*h);
    vx(i+1) = vx(i) + r1vx*h/2 + r2vx*h/2;
    vz(i+1) = vz(i) + r2vx*h/2 + r2vz*h/2;
    x(i+1) = x(i) + r1x*h/2 + r2x*h/2;
    z(i+1) = z(i) + r1z*h/2 + r2z*h/2;
    %parar quando atinge o solo
    if z(i) < 0
        break
    end
end 

figure(1)
plot(x(1:i),z(1:i),'y')         %primeiro o plot -> depois titulo, fundo, eixos, et
set(gca, 'color', 'k')
title('Trajetória do volante')
xlabel 'x (m)'
ylabel 'altura (m)'

figure(2)
plot(t(1:i),z(1:i),'y')         %primeiro o plot -> depois titulo, fundo, eixos, et
set(gca, 'color', 'k')
title('Evolução da altura do volante')
xlabel 'tempo (s)'
ylabel 'altura (m)'



%interpolacao para descobrir o ponto em que a função z(t) é 0
a = interp1(z(1:i),t(1:i),0)