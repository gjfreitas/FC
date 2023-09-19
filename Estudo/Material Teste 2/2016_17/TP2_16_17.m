%% 2o Teste Prático de Avaliação Discreta Física Computacional — 2016/2017
%% 1a)
clear all
close all
clc

L = 5;
h = 0.01;
x = 0:h:L;
N = length(x);

gama = -1.5;
eps = 2;

% Olhando para a equação temos uma matriz A:
% com -4e nas diagonais mas com 1 na primeira e ultima posição da diagonal
% 2e+h diagonalmente abaixo da diagonal execto (end,end-1) que é 0
% 2e-h diagonalmente acima da diagonal execto na posição (1,2) que é 0
% Temos também uma matriz b:
% (N,1) com 20 na primeira posição e 50 na ultima
% e preenchida por 2y(h^2) no resto da matriz

% Construção da Matriz A
A1 = diag(repmat(-4*eps,1,N));
A2 = diag(repmat(2*eps+h,1,N-1),-1);
A3 = diag(repmat(2*eps-h,1,N-1),1);
A = A1+A2+A3;

% É preciso fazer estes ajustes para garantir q T(1) = 20 e T(end) = 50
A(1,1) = 1;
A(1,2) = 0;
A(end,end) = 1;
A(end,end-1) = 0;

% Construção da Matriz b
b = repmat(2*h^2*gama,N,1);
b(1,1) = 20;
b(end,end) = 50;

Tmatrix = linsolve(A,b);
figure(1)
hold on
plot(x,Tmatrix,'b.-'), xlabel('x'), ylabel('T'), title('T em função de x')

% 1 b)
T_exact = 20 - gama.*x + (50 - 20 + gama*L)*(exp(x./eps)-1)/(exp(L/eps - 1));
plot(x,T_exact,'r.-')
legend('T obtido', 'T analitico')


%% 1 c)
clear all
close all
clc

L = 5;
h_k = [0.01 0.25 0.5];

for i = 1:length(h_k)
    h = h_k(i);
    x = 0:h:L;
    N = length(x);

    gama = -1.5;
    eps = 0.1;


    % Construção da Matriz A
    A1 = diag(repmat(-4*eps,1,N));
    A2 = diag(repmat(2*eps+h,1,N-1),-1);
    A3 = diag(repmat(2*eps-h,1,N-1),1);
    A = A1+A2+A3;
    A(1,1) = 1;
    A(1,2) = 0;
    A(end,end) = 1;
    A(end,end-1) = 0;

    % Construção da Matriz b
    b = repmat(2*h^2*gama,N,1);
    b(1,1) = 20;
    b(end,end) = 50;

    T = linsolve(A,b);
    
    if h == 0.01
        T1 = T;
        x1 = x;
    end
    if h == 0.25
        T2 = T;
        x2 = x;
    end
    if h == 0.5
        T3 = T;
    end
end
figure(1)
hold on
plot(x1,T1,'b.-',x2,T2,'r.-',x,T3,'m.-'), xlabel('x'), ylabel('T'), title('T em função de x (eps = 0.1)')
legend('T obtido h = 0.01', 'T obtido h = 0.25', 'T obtido h = 0.5')

% 1 b)
T_exact = 20 - gama.*x + (50 - 20 + gama*L)*(exp(x./eps)-1)/(exp(L/eps - 1));
figure(2)
plot(x,T_exact,'r.-',x,T3,'b.-'), xlabel('x'), ylabel('T'), title('T em função de x (eps = 0.1)')
legend('T obtido h = 0.5', 'T analitico h = 0.5')

%% 2 a)
clear all
close all
clc

p = 5:1:15;
N = 2.^p;
dw = 0.1;

wmax = (N.*dw)/2;
figure(1)
plot(p,log10(wmax),'r.-'), xlabel('p'), ylabel('log10(wmax)')

for i = 1:length(wmax)
    if wmax(i) >= 100
        break
    end
end

wm = wmax(i); % primeiro valor de Nyquist > 100
Nm = wm*2/dw; % valor de N correspondente
pm = log2(Nm); % valor de p corrspondente

dt = pi/wm;

t = 0:dt:((Nm-1)*dt);
w = 0:dw:((Nm-1)*dw);
y = cos(50*t) - 2*sin(90.1*t + pi/3);
figure(2)
plot(t,y), xlabel('t'), ylabel('y'), title('y = cos(50*t) - 2*sin(90.1*t + pi/3)')

%% 2 b)
clear all
close all
clc

% Usando dt = 0.025 e N = 2^13 => p = 13
dt = 0.025;
p = 13;
N = 2^p;
t = 0:dt:((N-1)*dt);
y = cos(50*t) - 2*sin(90.1*t + pi/3);

% Transformadas F(w)
Y = fft(y);
Y1 = fftshift(Y);

% Densidades Espectrais
DS = (dt*abs(Y)).^2;      % sem shift
DS1 = (dt*abs(Y1)).^2;    % com shift

% w com shift
% definir o eixo no espaço das frequências angulares para 0 ficar no centro
dw = 2*pi/(N*dt); % = delta omega
wmax = (N/2-1)*dw; % valor máximo no espaço das frequências
wmin = (-N/2)*dw;        % valor minimo no espaço das frequências
w(:,1) = wmin:dw:wmax;

plot(w,DS1,'r.'), title('y = cos(50t) - 2sin(90.1t + pi/3), com shift')

%% 2 c)
clear all
close all
clc

% Usando dt = 0.025 e N = 2^13 => p = 13
dt = 0.025;
p = 13;
N = 2^p;
t = 0:dt:((N-1)*dt);
y1 = cos(50*t) - 2*sin(90.1*t + pi/3);

for i = 1:length(t)
    Yt = dt * fftshift(fft(y1));      % transformada
    y2 = ifft(ifftshift(Yt)) / dt;   % transformada inversa
end
% Para obter a transformada de Fourier é preciso multiplicar por dt, pois o matlab calcula a transformada de fourier discreta
% slide 5 - Aula 7 (2020/21)

y = y1-y2;  
plot(t,y), xlabel('x'), ylabel('y1-y2'), title('Gráfico das diferenças entre y original e reconstruido')

mdiff = max(abs(y));
meandiff = mean(abs(y));
disp(['Maior diferença: ',num2str(mdiff)]);
disp(['Diferença média: ',num2str(meandiff)]);
disp('Como a maior diferença é um número bastante pequeno podemos que concluir que o metódo tem um erro quase nulo')