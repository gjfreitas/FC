%% Trabalho Prático 5
% Problemas de valor fronteira — Valores próprios
%% Problema 5.1: Método de shooting — determinação da frequência do primeiro modo normal de vibração 
clear all
close all
clc

u = 10^(-3); %kg/m
L = 1; %m
T = 10^3; %N

% calcular os primeiros 3 modos
n = [1 2 3];
wn = n*pi/L*sqrt(T/u);

w1 = wn(1); %rad/s
w2 = wn(2); %rad/s
w3 = wn(3); %rad/s

h = 0.01;
x = 0:h:L;
N = length(x);
y = nan(1,N);
y(1) = 0;
v = nan(1,N);
v(1) = 0.2; %m/s valor a escolha, diferente de 0

% i) Sabendo o w
for k = 1:N-1
   v(k+1) = v(k) + ((-w2^2 * u)/T * y(k))*h;
   y(k+1) = y(k) + v(k+1)* h ;
end
subplot(1,2,1)
plot(x,y), xlabel('x'),ylabel('y')
% O w escolhido  é a solução correta pois y(L) = y(0)

% ii) Usando um valor aleatório de w

wop = 3500;
for k = 1:N-1
   v(k+1) = v(k) + ((-wop^2*u)/T * y(k))*h;
   y(k+1) = y(k) + v(k+1)* h ;
end
subplot(1,2,2)
plot(x,y), xlabel('x'),ylabel('y')
% O w escolhido não é a solução correta pois y(L) é diferente do y(0)

%% Considerando W desconhecido

clear all
close all
clc

u = 10^(-3); %kg/m
L = 1; %m
T = 10^3; %N

h = 0.01;
x = 0:h:L;
N = length(x);

y = nan(N,1);
y(1) = 0;
v = nan(N,1);
v(1) = 2E-2;

%shooting
guess = [2000,3500]; % dois guesses inciais de w
result = [nan,nan];
B = 0; % o B corresponde ao valor pretendido de y(L), neste caso como queres y(L) = 0 assumimos B = 0

for k = 1:N-1
   v(k+1) = v(k) - ((guess(1).^2*u)/T * y(k))*h;
   y(k+1) = y(k) + v(k+1)* h ;
end

result(1) = y(end); % y do ultimo elemento

for k = 1:N-1
   v(k+1) = v(k) - ((guess(2).^2*u)/T * y(k))*h;
   y(k+1) = y(k) + v(k+1)* h ;
end

result(2) = y(end); % y do ultimo elemento

m = (result(2)-result(1))/(guess(2)-guess(1));
b = result(2) - m * guess(2);

guess(3) = guess(2)+ (B-result(2))/m;
for k = 1:N-1
   v(k+1) = v(k) - ((guess(3).^2*u)/T * y(k))*h;
   y(k+1) = y(k) + v(k+1)* h ;
end

result(3) = y(end); % y do ultimo elemento

while abs(guess(2)-guess(1)) > 1E-6 % O processo é repetido até a diferença dos guess não difira mais que uma determinada tolerância pré-estabelecida
    for k = 1:N-1
        v(k+1) = v(k) - ((guess(2).^2*u)/T * y(k))*h;
        y(k+1) = y(k) + v(k+1)* h ;
    end
    result(2) = y(end);
    m = (result(2)-result(1))/(guess(2)-guess(1));
    guess = [guess(2), guess(2)+ (B-result(2))/m];
    result(1) = result(2);
end
sol = guess(end-1);
disp(['omega: ', num2str(sol),' rad/s'])
figure(2)
plot(x,y,'.-'), xlabel('x'),ylabel('y')
grid

% um while de outra forma (comentar o while acima para este funcionar)
% i = 3;
% while abs(guess(i)-guess(i-1)) > 1E-6 % O processo é repetido até a diferença dos guess não difira mais que uma determinada tolerância pré-estabelecida
%     for k = 1:N-1
%         v(k+1) = v(k) - ((guess(i).^2*u)/T * y(k))*h;
%         y(k+1) = y(k) + v(k+1)* h ;
%     end
%     result(i) = y(end);
%     m = (result(i)-result(i-1))/(guess(i)-guess(i-1));
%     guess(i+1) = guess(i) + (B-result(i))/m;
%     
%     i = i+1;
% end
% sol = guess(end-1);
% disp(['omega: ', num2str(sol),' rad/s'])
% figure(1)
% plot(x,y,'.-'), xlabel('x'),ylabel('y')
% grid


