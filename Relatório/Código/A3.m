%% Trabalho Prático de Avaliação Contínua
% Gonçalo Freitas Nº Mec 98012
% Hugo Amaral Nº Mec 98224
% Turma PL7
%% A3)  Sobre-relaxação sucessiva (S.O.R).
%% Encontrar o valor de alfa optimo
clear all
close all
clc

L = 4;
h = 0.025;
N_max_iter = 1E6; % nº max de iterações
tol = 1E-6; % tolerância
delta_x = h;
delta_y = h;

y = -L/2:delta_y:L/2;
x = -L/2:delta_x:L/2;
N = length(x);
% [X,Y] = meshgrid(x,y);

% força atuadora
f = zeros(N,N);
for i = 1:length(x)
    for j = 1:length(y)
        if (x(i)^2 + y(j)^2) < (L^2)/9
            f(i,j) = 1;
        end
    end
end

% Condições de fronteira
z = zeros(N,N);
z(1,:) = 10;
z(end,:) = 12;
z(:,1) = 11 + 2/L .* x;
z(:,end) = 11 + 2/L .* x;


alfa_opt = 2 - (2*pi)/N; % alfa otimo teórico
alfa_vec = (alfa_opt-0.01):0.001:(alfa_opt+0.02); % vários valores de alfa perto do alfa optimo teórico

for ind_alfa = 1:length(alfa_vec)
    a = alfa_vec(ind_alfa);

    % Matriz z_old e z_new
    z_old = z;
    z_new = z_old;
    cond_fronteira = z_old == 0;
    
    for k = 1:N_max_iter
        for i = 2:(length(y)-1)
            for j = 2:(length(x)-1)
                if true(cond_fronteira(j,i))     
                    z_new(j,i) = (1-a)*z_old(j,i) + a/4*( z_new(j+1,i) + z_new(j-1,i) + z_new(j,i+1) + z_new(j,i-1) - (h^2)*f(j,i));
                end        
            end
        end

        dif = abs(z_new-z_old); % condição de convergência

        if max(dif) < tol
            break
        end

        z_old = z_new; % atualiza a estimativa da solução
    end
    iteracoes(ind_alfa) = k;
end

[u,v] = min(iteracoes);
disp(['Menor número de iterações : ',num2str(iteracoes(v))]);
disp(['Alfa para o menor número de iterações = ',num2str(alfa_vec(v)) ])
disp(['Alfa optimo teórico : ', num2str(alfa_opt)])
disp('---------------------------------------------------------')

% Agora que sabemos o alfa optimo podemos calcular a solução


% Matriz z_old e z_new
z_old = z;
z_new = z_old;

a = alfa_vec(v); % alfa otimo calculado

for k = 1:N_max_iter
    for i = 2:(length(y)-1)
        for j = 2:(length(x)-1)
            if true(cond_fronteira(j,i))     
                z_new(j,i) = (1-a)*z_old(j,i) + a/4*( z_new(j+1,i) + z_new(j-1,i) + z_new(j,i+1) + z_new(j,i-1) - (h^2)*f(j,i));
            end        
        end
    end
    
    dif = abs(z_new-z_old); % condição de convergência
    
    if max(dif) < tol
        break
    end
    
    z_old = z_new; % atualiza a estimativa da solução
end

figure
meshc(x,y,z_new), xlabel('x'), ylabel('y'), zlabel('z')
title('Método de Sobre-relaxação sucessiva')

figure
meshc(x,y,z_new), xlabel('x'), ylabel('y'), zlabel('z')
title('Método de Sobre-relaxação sucessiva')

disp(['Número de Iterações: ',num2str(k)])
