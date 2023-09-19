%% Trabalho Prático 8
%% Problema 8.1 Resolução de um sistema de equações lineares
%% 1 a) - Jacobi
clear all
close all
clc

N = 8;
% Construção da Matriz
A = diag(-ones([1 N]));
A(1,4) = sqrt(2)/2;
A(2,4) = sqrt(2)/2;
A(4,4) = -sqrt(2)/2;
A(4,6) = -1;
A(7,4) = -sqrt(2)/2;
A(1,5) = 1;
A(6,6) = 1;
A(3,7) = 1/2;
A(4,7) = 1/2;
A(7,7) = sqrt(3)/2;
A(8,7) = -sqrt(3)/2;
A(5,8) = 1;

b = zeros(N,1);
b(6) = 1E4;

x_old = ones(N,1);
x_new = x_old;

tol = 1E-7; % tolerância

kmax = 500; % nº max de iterações

for k = 1:kmax
    for i = 1:N
        x_new(i) = 1/A(i,i) * b(i); % ultimo termo
        soma = 0;
        for j = 1:N
            if j ~= i
                soma = soma + A(i,j)*x_old(j);
            end
        end
        x_new(i) = -1/A(i,i) * soma + x_new(i);
    end

    % tambem dava para fazer em dois ciclos for de 1 ate i-1 e de i+1 ate N, como na alinea b)
    
    dif = max(abs(x_new-x_old))/max(abs(x_new)); % condição de convergência
    
    if dif < tol
        break
    end
    
    x_old = x_new; % atualiza a estimativa da solução
end

disp(['Número de Iterações: ',num2str(k)])
disp(['Condição de Convergência: ',num2str(dif)])

%% 1 b) -  Gauss-Seidel
clear all
close all
clc

N = 8;
% Construção da Matriz
A = diag(-ones([1 N]));
A(1,4) = sqrt(2)/2;
A(2,4) = sqrt(2)/2;
A(4,4) = -sqrt(2)/2;
A(4,6) = -1;
A(7,4) = -sqrt(2)/2;
A(1,5) = 1;
A(6,6) = 1;
A(3,7) = 1/2;
A(4,7) = 1/2;
A(7,7) = sqrt(3)/2;
A(8,7) = -sqrt(3)/2;
A(5,8) = 1;

x_old = ones(N,1);
x_new = x_old;
b = zeros(N,1); b(6,1) = 1E4;

tol = 1E-7; % tolerância

kmax = 10000; % nº max de iterações

for k = 1:kmax
    for i = 1:N
        
        x_new(i) = 1/A(i,i) * b(i); % ultimo termo
        
        % j < i
        for j = 1:(i-1)
            x_new(i) = x_new(i) + (-1/A(i,i)) * A(i,j) * x_new(j);
        end
        
        % j > i
        for j = (i+1):N
            x_new(i) = x_new(i) + (-1/A(i,i)) * A(i,j) * x_new(j);
        end
    end
    
    dif = max(abs(x_new - x_old))/max(abs(x_new)); % condição de convergência
    
    if dif < tol
        break
    end
    
    x_old = x_new; % atualiza a estimativa da solução
end

disp(['Número de Iterações: ',num2str(k)])
disp(['Condição de Convergência: ',num2str(dif)])

%% 1 c) -  Sobre-relaxação sucessiva
clear all
close all
clc

alfa = 1.25;
N = 8;
% Construção da Matriz
A = diag(-ones([1 N]));
A(1,4) = sqrt(2)/2;
A(2,4) = sqrt(2)/2;
A(4,4) = -sqrt(2)/2;
A(4,6) = -1;
A(7,4) = -sqrt(2)/2;
A(1,5) = 1;
A(6,6) = 1;
A(3,7) = 1/2;
A(4,7) = 1/2;
A(7,7) = sqrt(3)/2;
A(8,7) = -sqrt(3)/2;
A(5,8) = 1;

x_old = ones(N,1);
x_new = x_old;
b = zeros(N,1); b(6,1) = 1E4;

tol = 1E-7; % tolerância

kmax = 10000; % nº max de iterações

for k = 1:kmax
    for i = 1:N
        
        x_new(i) = (1-alfa)*x_old(i); % ultimo termo
        soma = 0;
        for j = 1:N
            if j ~= i
                soma = soma + A(i,j)*x_new(j);
            end
        end
        x_new(i) = x_new(i) + alfa/A(i,i) * (-soma + b(i));
    end

    % tambem dava para fazer em dois ciclos for de 1 ate i-1 e de i+1 ate N, como na alinea b)
    
    dif = max(abs(x_new - x_old))/max(abs(x_new)); % condição de convergência
    
    if dif < tol
        break
    end
    
    x_old = x_new; % atualiza a estimativa da solução
end

disp(['Número de Iterações: ',num2str(k)])
disp(['Condição de Convergência: ',num2str(dif)])