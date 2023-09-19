%% 1o Teste PráticoFísica Computacional — 2015/2016 17 de março de 2016
clear all
close all
clc

% 1 a)
ro = 1.225; %kg/m^3
Cd = 0.508;
m = 57E-3; %kg
A = 3.5E-4; %m^2
z0 = 0.5; %m
v0 = 15; %m/s
g = 9.8; %m/s^2
h = 0.01; %s

t = 0:h:3; %s
N = length(t);
z = nan(1,N); z(1) = z0;
v = nan(1,N); v(1) = v0;

fv = @(t,V) - g - 0.5 * Cd * ro/m * A * abs(V)^2 * V;
fz = @(t,V) V;

hmax = z0;

for k = 1:N-1
    r1z = fz(t(k), v(k));
    r1v = fv(t(k), v(k));
    
    r2z = fz(t(k)+h/2, v(k)+r1v*h/2); 
    r2v = fv(t(k)+h/2, v(k)+r1v*h/2);
        
    r3z = fz(t(k)+h/2, v(k)+r2v*h/2); 
    r3v = fv(t(k)+h/2, v(k)+r2v*h/2);
    
    r4z = fz(t(k)+h, v(k)+r3v*h); 
    r4v = fv(t(k)+h, v(k)+r3v*h);
    
    z(k+1) = z(k) + 1/6 *(r1z + 2*r2z + 2*r3z +r4z)*h;
    v(k+1) = v(k) + 1/6 *(r1v + 2*r2v + 2*r3v +r4v)*h;

    %b)    
    if k > 1  
        if z(k) > z(k-1) && z(k) > z(k+1)
            hmax = z(k);
        end
    end
    
    %c)
    if z(k) < 0
        break
    end
end
figure(1)
plot(t,z), xlabel('t'), ylabel('z')
figure(2)
plot(t,v), xlabel('t'), ylabel('v')

% b)
% hmax = max(z);  %  outra maneira de achar o maximo
disp(['Altura máxima: ',num2str(hmax),' m'])

% c)
t_impacto = interp1([z(k-1),z(k)],[t(k-1),t(k)],0);
disp(['tempo de voo: ',num2str(t_impacto),' s'])


% zmax = [];
% for i=2:N-1
%     if z(i)>z(i-1) && z(i)>z(i+1)
%         aux = lagr(z(i-1:i+1), t(i-1:i+1));
%         zmax = [zmax, aux(1)]; 
%     end
% end 
% 
% max_z = zmax(1);
% disp(['Altura máxima: ', num2str(max_z), 'm'])
% 
% t_voo = interp1(z,t,0);
% disp(['Tempo de Voo: ', num2str(t_voo), 's'])

%% 1 d)
clear all
close all
clc

% 1 a)
ro0 = 1.225; %kg/m^3
Cd = 0.508;
m = 57E-3; %kg
A = 3.5E-4; %m^2
z0 = 7000; %m
v0 = 200; %m/s
g = 9.8; %m/s^2
h = 0.01; %s

t = 0:h:355; %s
N = length(t);
z = nan(1,N); z(1) = z0;
v = nan(1,N); v(1) = v0;

fv = @(t,Z,V) - g - 0.5 * Cd * ro0*exp(-Z/z0) * A/m * abs(V)^2 * V;
fz = @(t,V) V;


for k = 1:N-1
    
    
    r1z = fz(t(k), v(k));
    r1v = fv(t(k), z(k), v(k));
    
    r2z = fz(t(k)+h/2, v(k)+r1v*h/2); 
    r2v = fv(t(k)+h/2, z(k)+r1z*h/2, v(k)+r1v*h/2);
        
    r3z = fz(t(k)+h/2, v(k)+r2v*h/2); 
    r3v = fv(t(k)+h/2, z(k)+r1z*h/2, v(k)+r2v*h/2);
    
    r4z = fz(t(k)+h, v(k)+r3v*h); 
    r4v = fv(t(k)+h, z(k)+r1z*h, v(k)+r3v*h);
    
    z(k+1) = z(k) + 1/6 *(r1z + 2*r2z + 2*r3z +r4z)*h;
    v(k+1) = v(k) + 1/6 *(r1v + 2*r2v + 2*r3v +r4v)*h;
  
    %c)
    if z(k) < 0
        break
    end
end
figure(1)
plot(t,z), xlabel('t'), ylabel('z')
figure(2)
plot(t,v), xlabel('t'), ylabel('v')

% b)
hmax = max(z);  %  outra maneira de achar o maximo
disp(['Altura máxima: ',num2str(hmax),' m'])

% c)
t_impacto = interp1([z(k-1),z(k)],[t(k-1),t(k)],0);
disp(['tempo de voo: ',num2str(t_impacto),' s'])

% zmax = [];
% for i=2:N-1
%     if z(i)>z(i-1) && z(i)>z(i+1)
%         aux = lagr(z(i-1:i+1), t(i-1:i+1));
%         zmax = [zmax, aux(1)]; 
%     end
% end 
% 
% max_z = zmax(1);
% disp(['Valor máximo de "z": ', num2str(max_z), 'm'])
% 
% t_voo = interp1(z,t,0);
% disp(['Tempo de Voo: ', num2str(t_voo), 's'])


%% 2
% 2 a)
clear all
close all
clc

A = 8.3; %m/s^2
B = 21; %s^-1
z0 = 0.5; %m
v0 = 0; %m/s
h = 0.001; %s
t0 = 0; tf = 0.5; %s

t = t0:h:tf; %s
N = length(t);
z = nan(1,N); z(1) = z0;
v = nan(1,N); v(1) = v0;

A_CN = [0, 1 + B * h * 0.5; 1, -h/2];
b = [v0 * (1 - B * 0.5 * h) - 2 * A * h * 0.5; z0 + v0 * 0.5 * h];

for k = 1:(N-1)
    m_CN = linsolve(A_CN,b);
    z(k+1) = m_CN(1);
    v(k+1) = m_CN(2);
    
    b = [m_CN(2) * (1 - B * 0.5 *h) - 2 * A * h * 0.5; m_CN(1) + m_CN(2) * 0.5 * h];
end

figure(1)
plot(t,z), xlabel('t'), ylabel('z')
figure(2)
plot(t,v), xlabel('t'), ylabel('v')


%% 2 b)
clear all
close all
clc

A = 8.3; %m/s^2
B = 21; %s^-1
z0 = 0.5; %m
v0 = 0; %m/s
t0 = 0; tf = 0.5; %s
hs = [ 0.01 0.02 0.025 0.04 0.05];
erro = zeros(1, length(hs));

v_final = A/B * (exp(-B*tf) - 1);

for i = 1: length(hs)
    h = hs(i);
    t = t0:h:tf; %s
    N = length(t);
    z = nan(1,N); z(1) = z0;
    v = nan(1,N); v(1) = v0;
    
    A_CN = [0, 1 + B * h * 0.5; 1, -h/2];
    b = [v0 * (1 - B * 0.5 *h) - 2 * A * h * 0.5; z0 + v0 * 0.5 * h];

    for k = 1:(N-1)
        m_CN = linsolve(A_CN,b);
        z(k+1) = m_CN(1);
        v(k+1) = m_CN(2);

        b = [m_CN(2) * (1 - B * 0.5 *h) - 2 * A * h *0.5; m_CN(1) + m_CN(2) * 0.5 * h];
    end
    
    erro(i) = abs(v_final - v(end));
end
 
plot( log(hs), log(erro), '*b'), xlabel('log(h)'), ylabel('log(Erro)')
lsline

p = polyfit(log(hs), log(erro), 1);
ordem_do_metodo = p(1);         %ordem do metodo é o declive da reta
disp(['Ordem do metódo: ',num2str(round(ordem_do_metodo)),' '])

