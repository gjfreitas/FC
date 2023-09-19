%% Folha de Revisões
clear all
close all
clc

alfa = 0.2;
x_min = -20;
x_max = 20;
delta_x = 0.5;
x = x_min:delta_x:x_max;

z_min = 0;
z_max = 16;
delta_z = 0.1;

z = z_min:delta_z:z_max;

Nx = length(x);
Nz = length(z);

phi = zeros(Nx, Nz);

% Condições iniciais e fronteira
phi(:,1) = exp((-1/2).*x.^2);
phi(1,:) = 0;
phi(Nx,:) = 0;

figure(1)
plot(x,abs(phi(:,1)),'m.-'), xlabel('x(m)'), ylabel('abs(phi)')

eta = (1i*delta_z)/(4*(delta_x)^2);
csi = 2 * alfa * (delta_x)^2;

A2 = diag(-ones(Nx-3,1),1); % sobe um posição relativamente a diagonal
A3 = diag(-ones(Nx-3,1),-1);% desce uma posição relativamente a diagonal

A1 = diag(1/eta + 2 + csi.*x(2:(Nx-1)).^2);

A = A1+A2+A3;

% Crank-Nicolson
n = 1;
b = nan(Nx-2,1);

for j = 1:Nz-1
    for k = 1:Nx-2
        b(k,1) = phi(k,j) + ((1/eta)-2-csi*x(k+1).^2)*phi(k+1,j)+phi(k+2,j);
    end
    
    phi(2:Nx-1,j+1) = linsolve(A,b);
end

figure(2)
mesh(x,z,abs(phi)'), xlabel('x'), ylabel('z'), zlabel('abs(phi)')
figure(3)
contourf(x,z,abs(phi)'), xlabel('x'), ylabel('z')
h = colorbar; % barra de cores
set(get(h,'label'),'string','abs(phi)'); % label da barra de cores

% O método de Euler aplicado a esta equação seria sempre instável porque os 
% valores próprios da matriz são imaginários.