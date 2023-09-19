%% Trabalho Prático 2
%% Problema 2.1: Oscilador harmónico simples
clear all
close all 
clc

% e) Metódo de Euler-Cromer

K = 1; %N/m
m = 1; %kg
x0 =1; %m
v0 = 0; %m/s

t0 = 0; %s
tf = 50; %s
h = 0.01; %s

t = t0:h:tf; %s
N = length(t);

v = zeros(1,N);
v(1) = v0;

x = zeros(1,N);
x(1) = x0;

for k = 1:(N-1)
    v(k+1) = v(k) + ((-K/m)*x(k))*h;
    x(k+1) = x(k) + v(k+1)*h;  % unica mudança do Metódo de Euler-Cromer
end


count = 0;

for k = 2:(N-1) % começamos em 2 pois se começarmos em 1 não temos nenhum anterior para comparar (em matlab não posição 0)
    if x(k-1) <= x(k) && x(k) >= x(k+1)
        count = count + 1;
        
        aux = lagr(t(k-1:k+1),x(k-1:k+1));
        t_max(count) = aux(1);
        x_max(count) = aux(2);
    end
end

A = mean(x_max); %ampltiude
disp(['Amplitude média: ',num2str(A),' m'])

p = mean(diff(t_max)) ; %periodo
disp(['Período médio: ',num2str(p),' s'])

function lagr=lagr(xm,ym)
% determinacao de o ma'ximo de uma funcao discreta
%
% input: coordenadas de 3 pontos vizinhos de ordenadas maiores
%            matrizes xm e ym
% output: coordenadas do ponto ma'ximo (xmax,ymax)
%

xab=xm(1)-xm(2);
xac=xm(1)-xm(3);
xbc=xm(2)-xm(3);

a=ym(1)/(xab*xac);
b=-ym(2)/(xab*xbc);
c=ym(3)/(xac*xbc);

xml=(b+c)*xm(1)+(a+c)*xm(2)+(a+b)*xm(3);
xmax=0.5*xml/(a+b+c);

xta=xmax-xm(1);
xtb=xmax-xm(2);
xtc=xmax-xm(3);

ymax=a*xtb*xtc+b*xta*xtc+c*xta*xtb;

lagr(1)=xmax;
lagr(2)=ymax;
end