%% 2o Teste Prático de Avaliação Discreta Física Computacional — 2015/2016 26 de abril de 2016 — Teste A
%% 1 a)
clear all
close all
clc

beta = 18;
xi = 0;
xf = 7;
h = 0.01;
x = xi:h:xf;
N = length(x);

% Condições iniciais
Ti = 10E-4;
% Qi = 10E-4*v; %dT/dx = Q
T = nan(1,N); T(1) = Ti;
Q = nan(1,N);

B = 0;
guess = [1.9 2];
tol = 1E-6;
nshots = 100;
result = zeros(1,nshots);

for ishot = 1:nshots
    clear  Q
    Q(1) = Ti * guess(ishot);
    for k = 1:N-1
        Q(k+1) = Q(k) + (guess(ishot)*Q(k) - beta*exp(-1/T(k))*(1 + Q(k) - guess(ishot)*T(k)))*h;
        T(k+1) = T(k) + Q(k+1)*h; % Euler-Cromer
    end
    result(ishot) = Q(end);
    if ishot >= 2           % a partir da segunda iteração começamos a avaliar os results usando o metodo da secante
        m = (result(ishot)-result(ishot-1))/(guess(ishot)-guess(ishot-1));
        guess(ishot+1) = guess(ishot) + (B - result(ishot))/m;
        if abs(guess(ishot) - guess(ishot-1)) < tol
            break
       end      
    end
end
disp(['Velocidade: ',num2str(guess(ishot)),' m/s'])
plot(x,Q),xlabel('x'), ylabel('Derivada da Temperatura')
disp(['Derivada em x = 7: ',num2str(result(ishot))])


%% 1 b)
clear all
close all
clc

beta_k = 18:1:25;
xi = 0;
xf = 7;
h = 0.01;
x = xi:h:xf;
N = length(x);

% Condições iniciais
Ti = 10E-4;
% Qi = 10E-4*v; %dT/dx = Q
T = nan(length(beta_k),N); T(:,1) = Ti;
Q = nan(length(beta_k),N);

B = 0;
guess = [1.9 2.0];
tol = 1E-4;
nshots = 1001;
result = zeros(1,nshots);
v = nan(length(beta_k),2);
v(:,1) = beta_k(1,:); 

for j = 1:length(beta_k)
    beta = beta_k(j);
    clear guess
    guess = [1.9 2.0];
    for ishot = 1:nshots-1
        clear  Q
        Q(j,1) = Ti * guess(ishot);
        for k = 1:N-1
            Q(j,k+1) = Q(j,k) + (guess(ishot)*Q(j,k) - beta*exp(-1/T(j,k))*(1 + Q(j,k) - guess(ishot)*T(j,k)))*h;
            T(j,k+1) = T(j,k) + Q(j,k+1)*h; % Euler-Cromer
        end
        result(ishot) = Q(j,end);
        if ishot >= 2           % a partir da segunda iteração começamos a avaliar os results usando o metodo da secante
            m = (result(ishot)-result(ishot-1))/(guess(ishot)-guess(ishot-1));
            guess(ishot+1) = guess(ishot) + (B - result(ishot))/m;
            if abs(guess(ishot+1) - guess(ishot)) < tol
                v(j,2) = guess(ishot+1);
                break
           end      
        end
    end
end
figure(1)
plot(v(:,1),v(:,2)), xlabel('Beta'), ylabel('V')
% primeira coluna de v contem os valores de B
% segunda coluna de v contem o valor das velocidades

figure(2)
hold on
for i = 1:length(beta_k)
    plot(x,T(i,:))
end
xlabel('x')
ylabel('T')
legend('Beta = 18', 'Beta = 19', 'Beta = 20', 'Beta = 21', 'Beta = 23', 'Beta = 24', 'Beta = 25')
