%% 2º Teste Prático de Avaliação Discreta Física Computacional — 2018/2019
%% 1
clear all
close all
clc

% a) Não consegues ver freq menores que 1/(2*dt)
t0 = 0;
tf = 51.1;
N = 512;
dt = (tf-t0)/(N-1);
t = 0:dt:51.1;


% Função f(t)
y = 1.2*cos(5.0*t) + cos(5.1*t) + sin(5.5*t);

% 1 b)

% Transformadas F(w)
Y = fft(y);
YF = fftshift(Y);

DS = (dt*abs(YF)).^2;    % com shift


% w com shift
% definir o eixo no espaço das frequências angulares, >=0
dw = dt/(tf-t0);    % dw de modo a w ter o mesmo length que YF
wmax = 5.8;        % valor máximo no espaço das frequências
wmin = 4.8;        % valor minimo no espaço das frequências
w = wmin:dw:wmax;

plot(w,DS,'r.'), xlabel(' 4.8 <= w => 5.8'), ylabel('Densidade Espectral')

% 1 c)
YF = dt*fftshift(fft(y));
y2 = ifft(ifftshift(YF))/dt;   % transformada inversa
diffy = y-y2; 
figure(2)
plot(t,diffy,'b-',t,y,'k-',t,y2,'r-'), xlabel('t'), ylabel('y')
legend('Diferença', 'Y original', 'Y reconstruido')

mdiff = mean(abs(diffy));
maxdiff = max(abs(diffy));
disp(['Valor médio da diferença absoluta: ',num2str(mdiff)]);
disp(['Máximo da diferença absoluta: ',num2str(maxdiff)]);
disp('Como a maior diferença é um valor consideravel podemos concluir que existe alguma diferença entre a y original e y reconstruido')
