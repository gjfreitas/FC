%% Trabalho Prático 2
%% Problema 2.3 Oscilador quártico

clear all
close all
clc

% b)

m = 1; % kg
K = 1; % N/m
v0 = 0; % m/s
alfa = -0.1; % m^-2

t0 = 0; % s
h = 0.01; %s
tf = 50; %s

t = t0:h:tf;
N = length(t);

v = nan(1,N);
x = nan(1,N);
v(1) = v0;

x0_list = 0.1 : 0.1 : 2;


A_vec = nan(length(x0_list),1);
P_vec = nan(length(x0_list),1);

% Metódo de Euler-Cromer

for index_x0 = 1:length(x0_list)
    x(1) = x0_list(index_x0);
    for k = 1:(N-1)
        a = -K/m *(x(k) +2*alfa*x(k)^3);
        v(k+1) = v(k) + a*h;
        x(k+1) = x(k) + v(k+1)*h;%%Euler-cromer so faz com que aqui se meta mais 1
    end
    contador = 0;
    for k=2:N-1
        if x(k-1) <= x(k) && x(k) >= x(k+1)
            contador = contador+1;
            aux = lagr( t(k-1:k+1), x(k-1:k+1) );
            t_max(contador) = aux(1);
            x_max(contador) = aux(2);
        end
    end

    A = mean(x_max);
    P = mean(diff(t_max));
    clear t_max
    clear x_max
    A_vec(index_x0) = A;
    P_vec(index_x0) = P;
end
plot(A_vec,P_vec,'o-')

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