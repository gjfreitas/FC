clear all; close all; clc

% Pré-alocação
dz = 0.02;
dx = 0.1;
N = 1024;

x = (-(N-1)/2)*dx:dx:((N-1)/2)*dx; % x --> k

z0 = 0;
zf = 4;
z = z0:dz:zf;
Nz = length(z);

dk = 2*pi/(N*dx);
k = (-N/2)*dk:dk:((N/2)-1)*dk; % x --> k

q = zeros(Nz,N);
q(1,:) = sech(x);
Q = zeros(Nz,N);
Q(1,:) = fftshift(fft(q(1,:)));
qtexp0 = Q(1,:).*exp(1i.*k.^2*z(1)/2);

% ODE45
abstol = ones(1,N);
abstol = 10^(-9)*abstol;
options = odeset('RelTol',10^(-9),'AbsTol',abstol);

[z,qtexp] = ode45(@nonlinear, [z0:dz:zf], [qtexp0], options, N, k);

for i = 1:Nz
    q(i,:) = ifft(ifftshift(qtexp(i,:).*exp(1i.*k.^2.*z(i)./2)));
end


Amp = abs(q)*dx;
Int = abs(q).^2;

figure(1)
meshc(x,z,Int)
figure(2)
contourf(x,z,Int)
colorbar

