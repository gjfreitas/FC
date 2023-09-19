%% Folha de Revisões
clear all
close all
clc

load data.txt  %1024 pontos

Fs = 1000;            % Frequência de Amostragem                    
T = 1/Fs;             % Periodo de Amostragem       
N = length(data);     % Length of signal
t = (0:N-1)*T;        % Time vector

f = Fs*(0:(N/2))/N;

% centrar as frequências em x = 0 (espelhar as frequências)
f_sim = nan(N,1);
f_sim(1:(N/2)+1) = -f(end:-1:1);
f_sim((N/2):N) = f(1:1:end); 

zf = fft(data);
z = fftshift(zf);

dens_esp = (T*abs(z)).^2;


plot(f_sim,dens_esp), xlabel('f(Hz)'), ylabel('Densidade Espectral')

% b)
vmean = mean(dens_esp(dens_esp < 0.04)); % < 0.04 porque tudo o que tem densidade espectral menor que 0.04 é ruido
result = round(vmean*10^4)/10^4; % arredondar a 4 casas decimais
disp(['valor médio da densidade espetral do ruído = ',num2str(result)])


% % Outra forma de fazer b)
% for i = 1:length(dens_esp)
%     if dens_esp(i) < 0.04
%         ruido(i) = dens_esp(i);
%     end
% end
% 
% result = round(mean(ruido)*10^5)/10^5;% arredondar a 5 casas decimais
% disp(['valor médio da densidade espetral do ruído = ',num2str(result)])

% extra
figure(2)
plot(f_sim,abs(z)),xlabel('frequencia (Hz)'),ylabel('valor da transformada')

