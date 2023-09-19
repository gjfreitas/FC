function [t,ang,va] = rk2f(h,tf,angIni,va0,b,g,l,ome)

N = floor(tf/h);
va = zeros(N,1);
ang = zeros(N,1);
E = zeros(N,1);
ang(1) = angIni;
va(1) = va0;
t = 0;
fang = @(va) va;
fva = @(t,ang,va) -(g/l)*sin(ang)-b*va+10*sin(ome*t);


for k = 2:N;
    r1v = fva(t(k-1),ang(k-1),va(k-1));
    r1x = fang(va(k-1));
    r2v = fva(t(k-1)+h/2,ang(k-1)+(h/2)*r1x,va(k-1)+(h/2)*r1v);
    r2x = fang(va(k-1) + r1v*(h/2));
    va(k) = va(k-1) + r2v*h;
    ang(k) = ang(k-1) + r2x*h;
    t(k) = t(k-1) + (h/2);
end
end