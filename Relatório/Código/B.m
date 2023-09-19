%% Trabalho Prático de Avaliação Contínua
% Gonçalo Freitas Nº Mec 98012
% Hugo Amaral Nº Mec 98224
% Turma PL7
%% B1) o número de iterações é proporcional a M2 e B2) o tempo de cálculo é proporcional a M4
clear all
close all
clc

L = 4;
h_vec = [0.25 0.2 0.1 0.05 0.020 0.025];
HN = length(h_vec);
N_max_iter = 1E6; % nº max de iterações
tol = 1E-6; % tolerância

% pré alocação
N_vec = nan(1,HN);
k_vec = nan(1,HN); 
Time = nan(1,HN);

for u = 1:HN
    tic
    h = h_vec(u);

    delta_x = h;
    delta_y = h;

    y = -L/2:delta_y:L/2;
    x = -L/2:delta_x:L/2;
    N = length(x);
    N_vec(u) = N;

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

    % Matriz z_old e z_new

    z_old = z;
    z_new = z_old;


    cond_fronteira = z_old == 0;

    for k = 1:N_max_iter
        for i = 2:(length(y)-1)
            for j = 2:(length(x)-1)
                if true(cond_fronteira(j,i))     
                    z_new(j,i) = 1/4*( z_new(j+1,i) + z_new(j-1,i) + z_new(j,i+1) + z_new(j,i-1) - (h^2)*f(j,i));
                end        
            end
        end

        dif = abs(z_new-z_old); % condição de convergência

        if max(dif) < tol
            k_vec(u) = k;
            Time(u) = toc;
            break
        end
        
        z_old = z_new; % atualiza a estimativa da solução
    end
end

% B1)
figure(1)
plot(log(N_vec),log(k_vec)), xlabel('ln(M)'), ylabel('ln(Nº de iterações)')
p = polyfit(log(N_vec),log(k_vec),1);
Errop = (2-p(1))/2 * 100;
disp(['Declive ln(num iterações) em função de ln(M) : ',num2str(p(1)),'. Erro (considerando 2 como o valor esperado): ',num2str(Errop)])

% B2)
figure(2)
plot(log(N_vec), log(Time)), xlabel('ln(M)'), ylabel('ln(Tempo de cálculo) (s)')
q = polyfit(log(N_vec), log(Time),1);
Erroq = (4-q(1))/4 * 100;
disp(['Declive ln(tempo de cálculo) em função de ln(M) : ',num2str(q(1)),'. Erro (considerando 4 como o valor esperado): ',num2str(Erroq)])

