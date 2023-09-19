%% Problema 9.1 Poço de potencial infinito a uma dimensão
clear all
close all
clc

% S(x) = 0;
% g(x) = 2(E-V(x))
% y(x) = psi(x)

a = 1;
h = 0.001;

% En = n^2 * (pi^2)/(8*(a^2))
E1 = pi^2/(8*(a^2)); % n = 1
disp(['E1 : ',num2str(E1),' Ha'])

% Pré alocações
x = -a:h:a;
Nx = length(x);

psi = nan(1,Nx);
V = zeros(1,Nx); % V(x) = 0 dentro do potencial
g = nan(1,Nx);

% Condições fronteira
psi(1) = 0;
psi(2) = h*10^(-3); % previsão necessária pela fórmula do método de Numerov => psi(-a+h)

% guesses
guess(1) = 1.5;
guess(2) = 1.2;

% Goal do shooting
B = 0; % psi_f = 0

tol = 1E-12; % tolerância

for k = 2:Nx-1
    g(1:Nx) = 2*(guess(1)-V(1:Nx));
    psi(k+1) = (1 + h^2/12 * g(k+1))^(-1) * (-(1 + h^2/12 * g(k-1))*psi(k-1) + 2*(1 - 5 * h^2/12 * g(k))*psi(k));
end
result(1) = psi(end);

while abs(guess(2)-guess(1)) > tol
    g(1:Nx) = 2*(guess(1)-V(1:Nx));
    for k = 2:Nx-1
        psi(k+1) = (1 + h^2/12 * g(k+1))^(-1) * (-(1 + h^2/12 * g(k-1))*psi(k-1) + 2*(1 - 5 * h^2/12 * g(k))*psi(k));
    end
    result(2) = psi(end);
    m = (result(2)-result(1))/(guess(2)-guess(1));
    guess = [guess(2), guess(2) + (B-result(2))/m];
end
sol = guess(end-1);

C = trapz(x,abs(psi).^2);
psi_norm = psi/sqrt(C);
plot(x,psi_norm), xlabel('x'), ylabel('psi norm')

