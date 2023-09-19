%% Trabalho Prático 8
%% Problema 8.2 Resolução da equação de Laplace
%% 8.2 a) - Jacobi
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
% Exemplo: se h = 0.25, ind = [3,7], ou seja, entre os indices 3 e 7 a matriz potencial so tem 1's (ver Vi)

ind = [L/(2*h)+1 3*L/(2*h)+1];

% construir a matriz potencial (V_old)
V_old = zeros(N,N);
V_old(ind(1):ind(2),ind(1):ind(2)) = 1;
V_new = V_old;
Vi = V_old;

figure(1)
surf(X,Y,V_old)
title('Potencial inicial')

cond_fronteira = V_old == 0; % matriz booleana, 1 -> onde V_old == 0; 0 -> onde V_old == 1

for k = 1:kmax
    for i = 2:N-1
        for j = 2:N-1
            if true(cond_fronteira(i,j))
                V_new(j,i) = 1/4*( V_old(j+1,i) + V_old(j-1,i) + V_old(j,i+1) + V_old(j,i-1) );
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

figure
surf(X,Y,V_new)
title('Potencial Final')

disp(['Número de Iterações: ',num2str(k)])
disp(['Condição de Convergência: ',num2str(dif)])