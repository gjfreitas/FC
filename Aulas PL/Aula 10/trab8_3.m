%% Trabalho Prático 8
%% Problema 8.3 Cálculo do campo elétrico e da capacidade do condensador do problema 8.2
clear all
close all
clc

h = 0.025;  % usar h = 0.025 para dar um gráfico igual as soluções
tol = 1E-7; % tolerância
kmax = 2000; % nº max de iterações
L = 1;

x = -L:h:L;
y = -L:h:L;
[X,Y] = meshgrid(x,y);
N = length(x); % length(x) = length(y)


% a matriz tem valores zero, exceto no condutor interno
% ou seja x = -1/2 e +1/2; y = -1/2 e +1/2
% devemos identificar os indices correspondentes

ind = [L/(2*h)+1 3*L/(2*h)+1];

% construir a matriz potencial (V_old)
V_old = zeros(N,N);
V_old(ind(1):ind(2),ind(1):ind(2)) = 1;
V_new = V_old;
Vi = V_old;

% figure(1)
% surf(X,Y,V_old)
% title('Potencial inicial')

cond_fronteira = V_old == 0; % matriz booleana, 1 -> onde V_old == 0; 0 -> onde V_old == 1

for k = 1:kmax
    for i = 2:N-1
        for j = 2:N-1
            if true(cond_fronteira(i,j))
                V_new(i,j) = 1/4*(V_new(i+1,j) + V_new(i-1,j) + V_new(i,j+1) + V_new(i,j-1));
            end
        end
    end
    
    numerador = sqrt(sum(sum((V_new-V_old).^2)));
    denominadador = sqrt(sum(sum(V_new.^2)));
    dif = numerador/denominadador; % condição de convergência
    
    if dif < tol
        break
    end
    
    V_old = V_new; % atualiza a estimativa da soluçãp
end

figure(2)
surf(X,Y,V_new)
title('Potencial Final')

disp(['Número de Iterações: ',num2str(k)])
disp(['Condição de Convergência: ',num2str(dif)])

% Problema 8.3

[Ex, Ey] = gradient(V_new,h,h);
Ex = -Ex;
Ey = -Ey;
figure(3)
quiver(X,Y,Ex,Ey)

cap_1 = -trapz(x,Ex(:,1))*4; % capacidade pelo metódo 1 (carga)
display(['Capacidade pelo metódo 1: ',num2str(cap_1)])

E2 = Ex.^2 + Ey.^2;
cap_2 = mean(mean(E2))*4; % capacidade pelo metódo 2 (campo elétrico)
display(['Capacidade pelo metódo 2: ',num2str(cap_2)])
% Podem dar valores diferentes, tem a haver com o valor h

% Outra forma
% e0 = 1;
% Clambda_previsto = 2*pi/log(2); 
% deltaV = 1; % parede exterior tem V = 0 e parede interior tem V = 1;
% E = sqrt(Ex.^2 + Ey.^2);
% sigma = e0*E;
% parede1 = trapz(sigma(:,1),x);
% parede2 = trapz(sigma(:,N),x);
% parede3 = trapz(sigma(1,:),y);
% parede4 = trapz(sigma(N,:),y);
% Qlambda = abs(parede1 + parede2 + parede3 + parede4);
% Clambda = Qlambda/deltaV;