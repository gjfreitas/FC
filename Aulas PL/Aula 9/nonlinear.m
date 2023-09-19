function [ derivadas ] = nonlinear( z,qtexp,N,k )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

derivadas = zeros(N,1);
qt = qtexp./(exp(1i.*k.^2*z/2).');
q = ifft(ifftshift(qt));
V = 1i.*q.*q.*conj(q);
VT = fftshift(fft(V));
derivadas = VT.*(exp(1i.*k.^2*z/2).');
end
