%% Trabalho 1 - Parte 2
%% Problema 1.3 Movimento a três dimensões

clear all
close all 
clc

% f)
alpha = 40.6; % graus latitude N
W = 7.292E-5; % rad/s
w = [-W*cosd(alpha),0,W*sind(alpha)]; % rad/s
R_Terra = 6.37E6; % m
R = [0,0,R_Terra]; % m

% g) e h) considerando apenas o efeito da força centrífuga.
z0 = 200; % m
g = [0,0,-9.8]; %m/s^2
pos0 = [0,0,z0];
v0 = [0,0,0]; 
m = 1; %kg

t0 = 0; %s
h = 0.01; %s
tf = 10; %s
t = t0:h:tf; %s
ac = g - cross(w,cross(w,R));
 
N = length(t);
v = nan(N,3); % o mesmo que  v = zeros(2,N)
v(1,:) = v0;

pos = nan(N,3);
pos(1,:) = pos0;
 
 for k = 1:(N-1)
    
    if pos(k,3) < 0
        break
    end
    
    v(k+1,:) = v(k,:)+ac.*h;
    pos(k+1,:) = pos(k,:)+v(k,:).*h;
    
end

t_impacto = interp1([pos(k-1,3),pos(k,3)],[t(k-1),t(k)],0);
disp(['tempo de voo: ',num2str(t_impacto),' s'])
x_impacto = interp1([t(k-1),t(k)],[pos(k-1,1),pos(k,1)],t_impacto);
disp(['desvio em x: ',num2str(x_impacto),' m'])
y_impacto = interp1([t(k-1),t(k)],[pos(k-1,2),pos(k,2)],t_impacto);
disp(['desvio em y: ',num2str(y_impacto),' m'])

%% i) e j) considerando apenas o efeito da força de Coriolis.
clear all
close all 
clc


alpha = 40.6; % graus latitude N
W = 7.292E-5; % rad/s
w = [-W*cosd(alpha),0,W*sind(alpha)]; % rad/s
R_Terra = 6.37E6; % m
R = [0,0,R_Terra]; % m

z0 = 200; % m
g = [0,0,-9.8]; %m/s^2
pos0 = [0,0,z0];
v0 = [0,0,0]; 
m = 1; %kg

 t0 = 0; %s
 h = 0.01; %s
 tf = 10; %s
 t = t0:h:tf; %s
 
 N = length(t);
 v = nan(N,3); % o mesmo que  v = zeros(2,N)
 v(1,:) = v0;
 pos = nan(N,3);
 pos(1,:) = pos0;
 
 for k = 1:(N-1)
    
    if pos(k,3) < 0
        break
    end
    
    acorl = g - 2*cross(w,v(k,:));
    
    v(k+1,:) = v(k,:)+acorl.*h;
    pos(k+1,:) = pos(k,:)+v(k,:).*h;
    
end

t_impacto = interp1([pos(k-1,3),pos(k,3)],[t(k-1),t(k)],0);
disp(['tempo de voo: ',num2str(t_impacto),' s'])
x_impacto = interp1([t(k-1),t(k)],[pos(k-1,1),pos(k,1)],t_impacto);
disp(['desvio em x: ',num2str(x_impacto),' m'])
y_impacto = interp1([t(k-1),t(k)],[pos(k-1,2),pos(k,2)],t_impacto);
disp(['desvio em y: ',num2str(y_impacto),' m'])


%% k)
% - a força de Coriolis desvia para Este (e para Norte no hemisfério Sul);
% -a força centrífuga desvia para Sul (em ambos os hemisférios).