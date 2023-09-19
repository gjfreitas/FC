%% 1º Teste Prático de Avaliação Discreta Física Computacional — 2018/2019 29 de março de 2019
clear all
close all
clc

R = 10;
C = 1E-3;
eps = 5.0;
L = 0.10;

a = 1/(R*C) + R/L;
b = 2/(L*C);
c = eps/(L*C);

h = 0.000001;
tf = 10*R*C;
t = 0:h:tf;
N = length(t);
Vc = nan(1,N); Vc(1) = 0;
dVc = nan(1,N); dVc(1) = 0;

% 1 a)
for k = 1:N-1
    Vc(k+1) = Vc(k) + dVc(k)*h;
    dVc(k+1) = dVc(k) + (- a * dVc(k) - b * Vc(k) + c)*h;
    
    if abs(Vc(k)-Vc(k+1)) < 10^-6 && abs(eps/2 -Vc(k))< 10^-6
        break
    end

end
plot(t,Vc), xlabel('t(s)'), ylabel('Vc(V)')

for j=1:length(Vc)-1
    if Vc(j+1)<Vc(j) && Vc(j-1)<Vc(j)
        max=lagr(t(j-1:j+1),Vc(j-1:j+1));
    end
end
Vmax = max(2); % Vmax = max(V)
disp(['Máximo de Vc = ',num2str(Vmax),' V'])

%% 1 b) Euler implicito

clear all
close all
clc

R = 10;
C = 1E-3;
eps = 5.0;
L = 0.10;

a = 1/(R*C) + R/L;
b = 2/(L*C);
c = eps/(L*C);

dVc0 = 0;
Vc0 = 0;

h = 0.000001;
tf = 10*R*C;
t = 0:h:0.5;
N = length(t);
Vc = nan(1,N); Vc(1) = dVc0;
dVc = nan(1,N); dVc(1) = Vc0;

A = [ 1, -h; h*b, 1+a*h];
b = [ Vc0; dVc0 + c*h];

for k = 1:(N-1)
    z = linsolve(A,b);
    Vc(k+1) = z(1);
    dVc(k+1) = z(2);
    b = [z(1); z(2) + c*h];
    
    if abs(Vc(k)-Vc(k+1)) < 10^-6 && abs(eps/2 -Vc(k+1))< 10^-6
        break
    end
    
end 

plot(t,Vc,'.'), xlabel('t(s)'), ylabel('Vc(V)')

%% 1 c)
clear all
close all
clc

R = 10;
C = 1E-3;
eps = 5.0;
L = 0.10;

a = 1/(R*C) + R/L;
b = 2/(L*C);
c = eps/(L*C);

dVc0 = 0;
Vc0 = 0;

h=0.000001;
tf=10*R*C;
t = 0:h:0.5;
N = length(t);
Vc = nan(1,N); Vc(1) = dVc0;
dVc = nan(1,N); dVc(1) = Vc0;

fdv = @(Vc,dVc) - a * dVc - b * Vc + c;
fv = @(dVc) dVc;

for k=1:N-1
    r1dv = fdv( Vc(k),dVc(k) );
    r1v = fv(dVc(k));
    
    r2dv = fdv( Vc(k) + r1v * h/2 ,dVc(k) + r1dv*h/2 );
    r2v = fv(dVc(k) + r1dv*h/2 );
    
    Vc(k+1) = Vc(k) + r2v*h;
    dVc(k+1) = dVc(k) + r2dv*h;
    
    if abs(Vc(k)-Vc(k+1)) < 10^-6 && abs(eps/2 -Vc(k+1))< 10^-6
        break
    end
        
end

plot(t,Vc,'.'), xlabel('t(s)'), ylabel('Vc(V)')

%% 2 a) 2018/19
close all
clear all
clc

GmS=4*pi^2;x0=1.0167;vy0=8.2;
t0=0;tf=10;h=0.001;
t=t0:h:tf;
N=length(t);
x=zeros(2,N);x(1,1)=x0; 
v=zeros(2,N);v(2,1)=vy0;
theta=zeros(1,N); 
A=zeros(1,N);

for i=1:N-1
    r=norm(x(:,i));
    F=-GmS/r^3*x(:,i);
    v(:,i+1)=v(:,i)+h*F;
    x(:,i+1)=x(:,i)+h*v(:,i+1);
    theta(i+1)=mod(atan2(x(2,i),x(1,i)),2*pi);
    A(i)= (((norm(x(:,i+1))+norm(x(:,i)))/2)^2*(theta(i+1)-theta(i)))/2;
    if i~=1 && theta(i+1)<theta(i)
        theta(i+1)=theta(i+1)+2*pi;
        break
    end
end

T=interp1(theta(i-1:i+1),t(i-1:i+1),2*pi);
disp(['periodo: ',num2str(T),' anos'])
x=x(:,1:i);v=v(:,1:i);t=t(1:i);theta=theta(1:i);A=A(1:i);

figure(1)
plot(0,0,'oy','MarkerSize',10,'Linewidth',20);
hold on
plot(x(1,:),x(2,:),'.')
legend('Sol','órbita da Terra')
axis equal
% set(gca,'PlotBoxAspectRatio',[1 1 1])

%% 2 a) outra forma (2020/21)
clear all
close all
clc

h = 0.0001; %[ano]

t(:,1) = 0:h:1;
N = length(t);

x = nan(N,1);
y = nan(N,1);
r = nan(N,1);
ang = nan(N,1);
vx = nan(N,1);
vy = nan(N,1);


x(1) = 1.0167; % Unidades Astronomicas (AU)
y(1) = 0; % AU
r(1) = norm([x(1),y(1)]);
ang(1) = 0;
vx(1) = 0; % AU/ano
vy(1) = 8.2; % AU/ano

% Metódo de Euler-Cromer

for k = 1:(N-1)
    vx(k+1) = vx(k) - 4*pi^2* x(k) / (r(k)^3)*h;
    vy(k+1) = vy(k) - 4*pi^2* y(k) / (r(k)^3)*h;
    x(k+1) = x(k) + vx(k+1) * h;
    y(k+1) = y(k) + vy(k+1) * h;
    
    r(k+1) = norm([x(k+1),y(k+1)]);
    ang(k+1) = mod(atan2(y(k+1),x(k+1)),2*pi);
end

% a)
figure(1)
plot(x,y,'k.-'), xlabel('x'),ylabel('y')
set(gca,'PlotBoxAspectRatio',[1 1 1])

%b)

for k = 1 : N-1
    if ang(k+1) < ang(k)
        break
    end
end

disp(['periodo: ',num2str(t(k)),' anos'])

% c)

index_t_meio_periodo = floor(k/2);
area = nan(index_t_meio_periodo-1,1);

for k = 1:(index_t_meio_periodo-1)
    area(k) = ( ((r(k+1) + r(k)))^2)/2*(ang(k+1)-ang(k))/2;
end

figure(2)
plot(t(1:(index_t_meio_periodo-1)),area)
title('alinea c')
xlabel('t')
ylabel('Variação da Aréa')


