%% Problema 5.2: Método das diferenças finitas — determinação da frequência de vários modos normais de vibração
clear all
close all
clc

u = 10^(-3); %kg/m
L = 1; %m
T = 10^3; %N

h = 0.01;
x = 0:h:L;
N = length(x);

A1 = diag(repmat(-2,1,N-2)); % matriz com -2 na diagonal
A2 = diag(ones([1 N-3]),1); % sobe um posição relativamente a diagonal
A3 = diag(ones([1 N-3]),-1);% desce uma posição relativamente a diagonal


A = A1 + A2 + A3; % matriz pedida no enunciado

sol1 = eigs(A,3,'sm');

sol = sqrt(-sol1 * T / (u*h^2));
[vec,val]=eigs(A,3,'sm');

figure(1)
plot(plot(x(2:(end-1)),vec(:,1:3)),'.-')
