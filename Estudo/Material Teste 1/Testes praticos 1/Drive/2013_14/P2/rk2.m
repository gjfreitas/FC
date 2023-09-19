function [t,ang,va] = rk2(h,tf,angIni,va0,b,g,l)

N = floor(tf/h);
va = zeros(N,1);
ang = zeros(N,1);
E = zeros(N,1);
t = linspace(0,tf,N);
ang(1) = angIni;
va(1) = va0;

fang = @(va) va;
fva = @(ang,va) -(g/l)*sin(ang)-b*va;


for k = 2:N;
    r1v = fva(ang(k-1),va(k-1));
    r1x = fang(va(k-1));
    r2v = fva(ang(k-1)+(h/2)*r1x,va(k-1)+(h/2)*r1v);
    r2x = fang(va(k-1) + r1v*(h/2));
    va(k) = va(k-1) + r2v*h;
    ang(k) = ang(k-1) + r2x*h;
end
end