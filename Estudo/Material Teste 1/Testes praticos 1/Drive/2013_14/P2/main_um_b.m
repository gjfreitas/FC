clc, clear all, close all


% 1) b.

%
% Condi��es iniciais:
g = 9.8;
l = 1;
va0 = 0;
angIni = 0.2;
tf = 30;
h = 0.01;

% Varia��o de 'b' e c�lculo dos instantes correspondentes aos extremos
% cujo �ngulo � menor que (�ngulo inicial)/exp(1)) = 0.0736:
b = linspace(0.1,1,10);
for n = 1:length(b),
    [t,ang,va] = rk2(h,tf,angIni,va0,b(n),g,l);
    ic = 0;
    for q = 2:length(t),
        if va(q-1)*va(q)<0
        ic = ic + 1;
        ind(ic,n) = q;
        end
    end
    abs_ang = abs(ang);
    for a = 1:length(ind(:,n)),
        if abs_ang(ind(a,n)) < angIni/exp(1),
            tCONDICAO(a,n) = t(ind(a,n));
            break
        end
    end
end

% Instantes que correspondem � condi��o pedida.
% cada linha corresponde a um 'b' diferente, variando este entre 0,1 at� 1,
% com um passo de 0,1.
t_resposta = nonzeros(tCONDICAO)



