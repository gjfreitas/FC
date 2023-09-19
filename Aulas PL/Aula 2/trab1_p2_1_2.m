%% Trabalho 1 - Parte 2
%% Problema 1.2 Movimento a duas dimensões
 clear all
 close all
 clc
 
 % eixo dos zz aponta para cima e eixo dos xx para a direita
 % a) (d^2x)/(dt^2) = 0 e (d^2z)/(dt^2) = -g
 
 % b) e c) 
 
 theta = 37; % graus
 v0 = 50; %m/s
 v0x = v0*cosd(theta); %m/s
 v0z = v0*sind(theta); %m/s
 z0 = 55; %m
 m = 1; %kg
 g = 9.8; %m/s^2
 
 t0 = 0; %s
 h = 0.01; %s
 tf = 10; %s
 t = t0:h:tf; %s

N = length(t);
vx = nan(1,N); % o mesmo que  vx = zeros(1,N)
vx(1) = v0x;
vz = nan(1,N);
vz(1) = v0z;

x = nan(1,N);
x(1) = 0;
z = nan(1,N);
z(1) = z0;

for k = 1:(N-1) %N-1 para os comprimentos de v e t serem iguais
    
    if z(k) < 0
        break %penultimo valor de z positivo e ultimo negativo, isto faz-se para se poder usar interp1
    end
    
    vx(k+1) = vx(k);
    x(k+1) = x(k) + vx(k)*h;
    
    vz(k+1) = vz(k) - g*h;
    z(k+1) = z(k)+ vz(k)*h;
    
end

plot(t,z,'m-')
% plot(t,x,'k-') 
% plot(t,vz,'-b')
% plot(t,vx,'-r')

t_impacto = interp1([z(k-1),z(k)],[t(k-1),t(k)],0);
x_impacto = interp1([t(k-1),t(k)],[x(k-1),x(k)],t_impacto);
disp(['tempo de voo: ',num2str(t_impacto),' s'])
disp(['alcance: ',num2str(x_impacto),' m'])

%% b) e c) de outra forma, usando apenas uma matriz V em vez de VX e VY (forma vetorial em vez de forma escalar)
 clear all
 close all
 clc

 
 theta = 37; % graus
 v0 = 50; %m/s
 v0x = v0*cosd(theta); %m/s
 v0z = v0*sind(theta); %m/s
 z0 = 55; %m
 m = 1; %kg
 g = [0,0,-9.8]; %m/s^2
 
 t0 = 0; %s
 h = 0.01; %s
 tf = 10; %s
 t = t0:h:tf; %s

N = length(t);
v = nan(N,3); % o mesmo que  v = zeros(2,N)
v(1,:) = [v0x,0,v0z];

pos = nan(N,3);
pos(1,:) = [0,0,z0];

for k = 1:(N-1) %N-1 para os comprimentos de v e t serem iguais
    
    if pos(k,3) < 0
        break %penultimo valor de z positivo e ultimo negativo, isto faz-se para se poder usar interp1
    end
    
    v(k+1,:) = v(k,:)+g.*h;
    pos(k+1,:) = pos(k,:)+v(k,:).*h;
    
end

plot(t,pos(:,3),'m-')
t_impacto = interp1([pos(k-1,3),pos(k,3)],[t(k-1),t(k)],0);
disp(['tempo de voo: ',num2str(t_impacto),' s'])
x_impacto = interp1([t(k-1),t(k)],[pos(k-1,1),pos(k,1)],t_impacto);
disp(['alcance: ',num2str(x_impacto),' m'])


%% d) 
 clear all
 close all
 clc

 
 theta = 37; % graus
 v0 = 50; %m/s
 v0x = v0*cosd(theta); %m/s
 v0z = v0*sind(theta); %m/s
 z0 = 0; %m
 m = 1; %kg
 g = [0,0,-9.8]; %m/s^2
 
 t0 = 0; %s
 h = 0.01; %s
 tf = 10; %s
 t = t0:h:tf; %s

N = length(t);
v = nan(N,3); % o mesmo que  v = zeros(2,N)
v(1,:) = [v0x,0,v0z];

pos = nan(N,3);
pos(1,:) = [0,0,z0];

for k = 1:(N-1) %N-1 para os comprimentos de v e t serem iguais
    
    if pos(k,3) < 0
        break %penultimo valor de z positivo e ultimo negativo, isto faz-se para se poder usar interp1
    end
    
    v(k+1,:) = v(k,:)+g.*h;
    pos(k+1,:) = pos(k,:)+v(k,:).*h;
    
end


t_impacto = interp1([pos(k-1,3),pos(k,3)],[t(k-1),t(k)],0);
disp(['tempo de voo: ',num2str(t_impacto),' s'])
x_impacto = interp1([t(k-1),t(k)],[pos(k-1,1),pos(k,1)],t_impacto);
disp(['alcance: ',num2str(x_impacto),' m'])
r = ((v0)^2/norm(g))*sind(2*theta);
disp(['alcance exacto: ',num2str(r),' m'])
erro = abs(x_impacto-r)/r *100;
disp(['erro: ',num2str(erro),' %'])