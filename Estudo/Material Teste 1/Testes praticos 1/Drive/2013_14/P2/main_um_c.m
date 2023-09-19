clear all, close all, clc

% Condições iniciais:
% A massa é 1.
% Omega = 1.
g = 9.8;
l = 1;
b = 1;
va0 = 0;
angIni = 0.2;
tf = 120;
h = 0.01;
ome = 1;


[t,ang,va] = rk2f(h,tf,angIni,va0,b,g,l,ome);

% Indexação dos pontos onde a velocidade troca de sentido e interpolação:
ic = 0;
for n = floor(length(t)/2):length(t),
    if va(n-1)*va(n)<0
        ic = ic+1;
        ind(ic) = n;
    end
end
for n = 1:ic,
    tI1(n,:) = t(ind(n)-1:ind(n)+1);
    angI1(n,:) = ang(ind(n)-1:ind(n)+1);
    vaI1(n,:) = va(ind(n)-1:ind(n)+1);
end
for n = 1:ic,
    tN1(n) = interp1(vaI1(n,:),tI1(n,:),0);
    angN1(n) = interp1(vaI1(n,:),angI1(n,:),0);
end

% Cálculo do período através da média do tempo entre picos.
Periodo1 = mean(diff(tN1))*2

% Cálculo da amplitude através da média da posição entre extremos positivos
% e negativos.
for j = 1:length(angN1)/2,
    angSN1(j) = abs(angN1(2*j)-abs(angN1(2*j-1)));
end
Amplitude1 = mean(angSN1)/2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
% Omega = 2;
% Condições iniciais:
% A massa é 1.
g = 9.8;
l = 1;
b = 1;
va0 = 0;
angIni = 0.2;
tf = 120;
h = 0.01;
ome = 2;


[t,ang,va] = rk2f(h,tf,angIni,va0,b,g,l,ome);

% Indexação dos pontos onde a velocidade troca de sentido e interpolação:
ic = 0;
for n = floor(length(t)/2):length(t),
    if va(n-1)*va(n)<0
        ic = ic+1;
        ind(ic) = n;
    end
end
for n = 1:ic,
    tI2(n,:) = t(ind(n)-1:ind(n)+1);
    angI2(n,:) = ang(ind(n)-1:ind(n)+1);
    vaI2(n,:) = va(ind(n)-1:ind(n)+1);
end
for n = 1:ic,
    tN2(n) = interp1(vaI2(n,:),tI2(n,:),0);
    angN2(n) = interp1(vaI2(n,:),angI2(n,:),0);
end

% Cálculo do período através da média do tempo entre picos.
Periodo2 = mean(diff(tN2))*2

% Cálculo da amplitude através da média da posição entre extremos positivos
% e negativos.
for j = 1:length(angN2)/2,
    angSN2(j) = abs(angN2(2*j)-abs(angN2(2*j-1)));
end
Amplitude2 = mean(angSN2)/2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
% Omega = 3.
% Condições iniciais:
% A massa é 1.
g = 9.8;
l = 1;
b = 1;
va0 = 0;
angIni = 0.2;
tf = 120;
h = 0.01;
ome = 3;


[t,ang,va] = rk2f(h,tf,angIni,va0,b,g,l,ome);
% Indexação dos pontos onde a velocidade troca de sentido e interpolação:
ic = 0;
for n = floor(length(t)/2):length(t),
    if va(n-1)*va(n)<0
        ic = ic+1;
        ind(ic) = n;
    end
end
for n = 1:ic,
    tI3(n,:) = t(ind(n)-1:ind(n)+1);
    angI3(n,:) = ang(ind(n)-1:ind(n)+1);
    vaI3(n,:) = va(ind(n)-1:ind(n)+1);
end
for n = 1:ic,
    tN3(n) = interp1(vaI3(n,:),tI3(n,:),0);
    angN3(n) = interp1(vaI3(n,:),angI3(n,:),0);
end

% Cálculo do período através da média do tempo entre picos.
Periodo3 = mean(diff(tN3))*2

% Cálculo da amplitude através da média da posição entre extremos positivos
% e negativos.
for j = 1:length(angN3)/2,
    angSN3(j) = abs(angN3(2*j)-abs(angN3(2*j-1)));
end
Amplitude3 = mean(angSN3)/2

