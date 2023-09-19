%% 1o Teste Prático Física Computacional — 2011/2012 20 de Março de 2012
% Turma P1

%% 1 96P/Machholz

clear all
close all
clc

x0 = -0.1246177; %AU
v0y = -24.9099; %AU/ano
y0 = 0; %AU
v0x = 0; %AU/ano
r0=[x0 y0];
v0=[v0x v0y];

% a = G*ms = 4pi^2

a=4*pi^2;

h=0.0001;
t=0:h:10;
N=length(t);
x = nan(N,1); x(1) = x0;
y = nan(N,1); y(1) = y0;
r = nan(N,1); r(1) = norm([x(1),y(1)]);
ang = nan(N,1); ang(1) = 0;
vx = nan(N,1); vx(1) = v0x;
vy = nan(N,1); vy(1) = v0y;
r(1) = norm([x(1),y(1)]);
v(1,:)=v0;
t(1)=0;



for k=1:N-1
    vx(k+1) = vx(k) - 4*pi^2* x(k) / (r(k)^3)*h;
    vy(k+1) = vy(k) - 4*pi^2* y(k) / (r(k)^3)*h;
    x(k+1) = x(k) + vx(k+1) * h;
    y(k+1) = y(k) + vy(k+1) * h;
    
    r(k+1) = norm([x(k+1),y(k+1)]);
    ang(k+1) = mod(atan2(y(k+1),x(k+1)),2*pi);
        
end


plot(x,y,'b.'), xlabel('x (UA)'), ylabel('y (UA)')
grid on


%  b)
for i = 1 : N-1
    if ang(i+1) < ang(i)
        break
    end
end

disp(['periodo: ',num2str(t(i)),' anos'])
disp(['distância máxima: ',num2str(max(r)),' AU'])

% c)

index_t_meio_periodo = floor(i/2);
area = nan(index_t_meio_periodo-1,1);

for j = 1:(index_t_meio_periodo-1)
    area(j) = ( ((r(j+1) + r(j)))^2)/2*(ang(j+1)-ang(j))/2;
end

disp(['Area : ',num2str(area(j)),' AU^2'])

% d)
periodo = t(i);
t2=0:h:periodo/2;
N2=length(t2);
r2=zeros(N2,2);
v=zeros(N2,2);
r2(1,:)=r0;
v2(1,:)=v0;

% a = G * ms = 4 * pi^2
% Ec/Ep = (1/2 * mc * v^2)/(-a*m/r) = (-v^2 * r)/ 2a

Er=zeros(N2,1);
Er(1)=-((v2(1))^2*norm(r2(1,:)))/(2*a);

for k=1:N2-1
    
    v2(k+1,:)=v2(k,:)-(a/(norm(r2(k,:)))^3)*r2(k,:)*h;
    r2(k+1,:)=r2(k,:)+v2(k+1,:)*h;
    Er(k+1)=-((v2(k+1))^2*norm(r2(k+1,:)))/(2*a);
end

figure
plot(t2,Er,'b.'), xlabel('t (ano)'), ylabel('Ec/Ep (J)')
grid on

%% 2 Hale-Bopp
clear all
close all
clc

% 1 a)
x0 = -0.9141335; %AU
v0y = -9.2822965; %AU/ano
y0 = 0; %AU
v0x = 0; %AU/ano
r0=[x0 y0];
v0=[v0x v0y];

%a=G*ms
a=4*pi^2;

N=4000000;
x = nan(N,1); x(1) = x0;
y = nan(N,1); y(1) = y0;
r = nan(N,1); r(1) = norm([x(1),y(1)]);
ang = nan(N,1); ang(1) = 0;
vx = nan(N,1); vx(1) = v0x;
vy = nan(N,1); vy(1) = v0y;
r(1) = norm([x(1),y(1)]);
v(1,:)=v0;
t(1)=0;



for k=1:N-1
    
    if norm(r(k,:))<=70
        h=0.00005;
    else 
        h=0.0025;
    end
    
    t(k+1)=t(k)+h;
    
    vx(k+1) = vx(k) - 4*pi^2* x(k) / (r(k)^3)*h;
    vy(k+1) = vy(k) - 4*pi^2* y(k) / (r(k)^3)*h;
    x(k+1) = x(k) + vx(k+1) * h;
    y(k+1) = y(k) + vy(k+1) * h;
    
    r(k+1) = norm([x(k+1),y(k+1)]);
    ang(k+1) = mod(atan2(y(k+1),x(k+1)),2*pi);
        
end


plot(x,y,'b.'), xlabel('x (UA)'), ylabel('y (UA)')
grid on

% b)
for k = 1 : N-1
    if ang(k+1) < ang(k)
        break
    end
end

periodo = t(k);
disp(['periodo: ',num2str(periodo),' anos'])