%% 1o Teste Prático Física Computacional — 2013/2014 11 março de 2014
%Turma P1
clear all
close all
clc

g = 9.8; %m/s^2
r = 0.03; %m
x0 = 2.5; %m
V0 = (pi/3)*(15.625); %m^3
A0 = 15.625/x0; %m^2

% a)
R = sqrt(A0/pi);

theta = atan(R/x0); %rad
disp(['Ângulo do tanque cónico: ',num2str(theta),' rads'])

h = 2;

tf = 8*60; %seg

% b)
t = 0:h:tf;
N = length(t);
x = zeros(N,1);
x(1)= x0;

% A(x)=pi*(tan(theta))^2*x^2

f=@(x) -0.6*pi*r^2*(sqrt(2*g*x)/(pi*(tan(theta))^2*x^2));

for k=1:N-1
    
    r1 = f(x(k));
    r2 = f(x(k)+r1*h/2);
    r3 = f(x(k)+r2*h/2);
    r4 = f(x(k)+r3*h/2);
    
    x(k+1) = x(k)+(1/6)*(r1+2*r2+2*r3+r4)*h;
    
    if x(k+1)<0
        break
    end
    
end

plot(t,x,'r.'), xlabel('t'), ylabel('x')
grid on

%% c)
clear all
close all
clc

g = 9.8; %m/s^2
r = 0.03; %m
x0 = 2.5; %m
V0 = (pi/3)*(15.625); %m^3
A0 = 15.625/x0; %m^2


R = sqrt(A0/pi);

theta = atan(R/x0); %rad
disp(['Ângulo do tanque cónico: ',num2str(theta),' rads'])

hs = 2:0.5:4 ;

tf = 8*60; %seg

for i = 1: length(hs)
    h = hs(i);
    t = 0:h:tf;
    N = length(t);
    x = zeros(N,1);
    x(1)= x0;

    % A(x)=pi*(tan(theta))^2*x^2

    f=@(x) -0.6*pi*r^2*(sqrt(2*g*x)/(pi*(tan(theta))^2*x^2));
    for k=1:N-1

        r1 = f(x(k));
        r2 = f(x(k)+r1*h/2);
        r3 = f(x(k)+r2*h/2);
        r4 = f(x(k)+r3*h/2);

        x(k+1) = x(k)+(1/6)*(r1+2*r2+2*r3+r4)*h;

        if x(k+1)<0
            break
        end

    end
end
plot(t,x,'r.'), xlabel('t'), ylabel('x')
grid on

disp(['Valor esperado: ',num2str(x(end)),' m'])


%% d)

clear all
close all
clc

g = 9.8; %m/s^2
r = 0.03; %m
x0 = 2.5; %m
V0 = (pi/3)*(15.625); %m^3
A0 = 15.625/x0; %m^2


R = sqrt(A0/pi);

theta = atan(R/x0); %rad
disp(['Ângulo do tanque cónico: ',num2str(theta),' rads'])

hs = 2 ;

tf = 8.7*60; %seg

for i = 1: length(hs)
    h = hs(i);
    t = 0:h:tf;
    N = length(t);
    x = zeros(N,1);
    x(1)= x0;

    %A(x)=pi*(tan(theta))^2*x^2

    f=@(x) -0.6*pi*r^2*(sqrt(2*g*x)/(pi*(tan(theta))^2*x^2));
    for k=1:N-1

        r1 = f(x(k));
        r2 = f(x(k)+r1*h/2);
        r3 = f(x(k)+r2*h/2);
        r4 = f(x(k)+r3*h/2);

        x(k+1) = x(k)+(1/6)*(r1+2*r2+2*r3+r4)*h;

        if x(k+1) == 0
            ind = k;
            break
        end
        
        if x(k+1) < 0
            break
        end

    end
end
plot(t,x,'r.'), xlabel('t'), ylabel('x')
grid on

disp(['Valor esperado: ',num2str(x(end)),' m'])