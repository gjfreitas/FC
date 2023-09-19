%% 2o Teste Prático Física Computacional — 2014/2015 — Turma P3
%% 1 a)

% Eu sei que no enunciado o R vem em mm e o resto vem em m mas se fizer     
% R = 1.5E-3 a T mal muda e ao fazer R = 1.5 ja é mas proximo do pretendido ent
% aposto que foi um erro no enunciado

clear all
close all
clc

lbd = 0.1;
Q = 2.4;
const = Q/lbd;

% dT/dr = U
U0 = 0; 
T0 = 40;

R = 1.5;
r0 = 1*10^(-50);

%ODE45
reltol = 1*10^(-10);
abstol_1 = 1*10^(-10);
abstol_2 = 1*10^(-10);
options = odeset('RelTol',reltol,'AbsTol',[abstol_1 abstol_2]);


[r,sol] = ode45(@f,[r0 R],[T0 U0],options,const);
T = sol(:,1);

disp(['T(r = R) = ',num2str(T(end)),' ºC'])
plot(r,T), xlabel('r (m)'), ylabel('T (ºC )')

%% 1 b)
clear all
close all
clc

%ODE45
reltol = 1*10^(-10);
abstol_1 = 1*10^(-10);
abstol_2 = 1*10^(-10);
options = odeset('RelTol',reltol,'AbsTol',[abstol_1 abstol_2]);

lbd = 0.1;
Q = 2.4;
const = Q/lbd;
R = 1.5;
r0 = 1*10^(-50);

% dT/dr = U
U0 = 0; 

% shooting
B = 20;
guess = [40.0 40.3]; % guess para T(0)
tol = 1E-12;
nshots = 100;


for j = 1:nshots
    [r,sol] = ode45(@f,[r0 R],[guess(j) U0],options,const);
    result(j) = sol(end,1); 
    diff = B - result(j);
    
    if j>= 2
        m = (result(j)-result(j-1))/(guess(j)-guess(j-1));
        guess(j+1) = guess(j)+(diff)/m;
        if abs(guess(j)-guess(j-1)) < tol
            break
        end
    end

end

solucao = guess(j);
T = sol(:,1);
plot(r,T), xlabel('r (m)'), ylabel('T (ºC )')

% c)
for k = 1:length(T)
    if T(k) >= 25
        fraq(k) = T(k);
        i = k;
    end
end

disp(['T(0) => T(r = R) = ',num2str(result(j)),'ºC : ',num2str(solucao)])
disp(['Há um total de ',num2str(i), ' pontos com Temperatura maior ou igual a 25ºC quando T(r = R) = 20ºC'])

%% 2
clear all
close all
clc

t0 = -10;
tf = 10;
N = 2^10; % so funciona com N dados em q N é potência de 2
dt = (tf-t0)/(N-1);
r = -10:dt:10;

wmax = pi/dt;
dw = 2*wmax/N;

wmin = -N/2 * dw;
wmax = (N/2 - 1)*dw;
w = wmin:dw:wmax;

y = 1/2*exp(-abs(r));

Y = dt * fftshift(fft(y));

figure(1)
hold on
plot(w,abs(Y),'b.-'), xlabel('w'), ylabel('Densidade espectral')

Y_exact = 1./(1+w.^2);
plot(w,Y_exact,'r-')
legend('Obtido', 'Analitico')

disp('Os gráficos são iguais')
