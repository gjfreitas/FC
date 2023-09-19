%% Trabalho Prático 7 - Transformada de Fourier discreta e sua aplicação na resolução de equações diferenciais
%% Problema 7.2: Resolução da equação paraxial
clear all
close all
clc

dx = 0.05; % passo de x
N = 2^10; % numeros de pontos em x
x = (-(N-1)/2)*dx:dx:((N-1)/2)*dx;

z0 = 0;
dz = 0.02;
zf = 4;
z = z0:dz:zf;
Nz = length(z);

q = zeros(Nz,N); % pré-alocação
q(1,:) = exp(-x.^2./2); % a) forma do feixe em z = 0
% q(1,:) = sech(x); % b) forma do feixe em z = 0 -> Secante hiperbólica

dk = 2*pi/(N*dx);
kmax = (N/2-1)*dk;
kmin = (-N/2)*dk;
k = kmin:dk:kmax;

% Transformada de Fourier para z = 0; centrada
Qt0 = fftshift(fft(q(1,:)));

Qt = zeros(Nz,N); % pré-alocação

for i = 1:Nz
    Qt = Qt0.*exp(-1i.*k.^2.*z(i)/2);
    q(i,:) = ifft(ifftshift(Qt));
end

figure(1)
perfil = abs(q).^2;
mesh(x,z,perfil), xlabel('x'), ylabel('z'), zlabel('Intensidade |q|^2')
figure(2)
contourf(x,z,perfil), xlabel('x'), ylabel('z')  % T' transposta de T para dar igual as soluções
h = colorbar; % barra de cores
set(get(h,'label'),'string','Intensidade |q|^2'); % label da barra de cores
