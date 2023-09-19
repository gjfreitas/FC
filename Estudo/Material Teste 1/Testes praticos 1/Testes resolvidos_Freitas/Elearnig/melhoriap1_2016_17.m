%% 1o Teste Prático — Melhoria Física Computacional — 2016/2017 23 de junho de 2017
clear all
close all
clc

% 1 a)
x0 = -0.9141335; %AU
v0y = -9.2822965; %AU/ano
y0 = 0; %AU
v0x = 0; %AU/ano
r0=[x0 y0];
v0=[v0x v0y];

%a=G*ms
a=4*pi^2;

N=4000000;
x = nan(N,1); x(1) = x0;
y = nan(N,1); y(1) = y0;
r = nan(N,1); r(1) = norm([x(1),y(1)]);
ang = nan(N,1); ang(1) = 0;
vx = nan(N,1); vx(1) = v0x;
vy = nan(N,1); vy(1) = v0y;
r(1) = norm([x(1),y(1)]);
t = nan(N,1);
t(1)=0;



for k=1:N-1
    
    if norm(r(k,:))<=70
        h=0.00005;
    else 
        h=0.0025;
    end
    
    t(k+1)=t(k)+h;
    
    vx(k+1) = vx(k) - 4*pi^2* x(k) / (r(k)^3)*h;
    vy(k+1) = vy(k) - 4*pi^2* y(k) / (r(k)^3)*h;
    x(k+1) = x(k) + vx(k+1) * h;
    y(k+1) = y(k) + vy(k+1) * h;
    
    r(k+1) = norm([x(k+1),y(k+1)]);
    ang(k+1) = mod(atan2(y(k+1),x(k+1)),2*pi);
        
end


plot(x,y,'b.'), xlabel('x(UA)'), ylabel('y(UA)')
grid on
axis([-100 400 -400 400])

% b)
for k = 1 : N-1
    if ang(k+1) < ang(k)
        break
    end
end

periodo = t(k);
disp(['periodo: ',num2str(periodo),' anos'])

% c)
index_t_meio_periodo = floor(k/2);
area = nan(index_t_meio_periodo-1,1);

for k = 1:(index_t_meio_periodo-1)
    area(k) = ( ((r(k+1) + r(k)))^2)/2*(ang(k+1)-ang(k))/2;
end
disp(['Area : ',num2str(area(k)),' AU^2'])

% d)
% a = G * ms = 4 * pi^2
% Ec/Ep = (1/2 * mc * v^2)/(-a*m/r) = (-v^2 * r)/ 2a

t2 = 0:h:periodo/2;
N2 = length(t2);
r2 = zeros(N2,2);
v2 = zeros(N2,2);
r2(1,:) = r0;
v2(1,:) = v0;
Er = zeros(N2,1);
Er(1) = -((v2(1))^2*norm(r2(1,:)))/(2*a);

for k=1:N2-1
    
    v2(k+1,:)=v2(k,:)-(a/(norm(r2(k,:)))^3)*r2(k,:)*h;
    r2(k+1,:)=r2(k,:)+v2(k+1,:)*h;
    Er(k+1)=-((v2(k+1))^2*norm(r2(k+1,:)))/(2*a);
end

figure
plot(t2,Er,'b.'), xlabel('t(s)'), ylabel('Ec/Ep (J)')
grid on

%% 2 a)
clear all
close all
clc

K = 1; %N/m
m = 2; %Kg
x0 = 1; %m
v0 = 0; %m/s

h = 0.01;
t = 0:h:100;
N = length(t);

x = nan(1,N); x(1) = x0;
v = nan(1,N); v(1) = v0;

fv = @(X) -K/m * X;
fx = @(V) V;

for k=1:N-1
    
    r1x = fx(v(k));
    r1v = fv(x(k));
    
    r2x = fx(v(k)+r1v*h/2);
    r2v = fv(x(k)+r1x*h/2);
    
    r3x = fx(v(k)+(2*r2v-r1v)*h);
    r3v = fv(x(k)+(2*r2x-r1x)*h);
    
    x(k+1) = x(k) +1/6 * (r1x + 4*r2x + r3x) * h;
    v(k+1) = v(k) +1/6 * (r1v + 4*r2v + r3v) * h;
end
figure(1)
plot(t,x), xlabel('t(s)'), ylabel('x(m)')

Ec = 0.5*m*v.^2;
Ep = 0.5*K*x.^2;
Em = Ec + Ep;
figure(2)
plot(t,Em), xlabel('t(s)'), ylabel('Em(J)')

Em_i = 0.5 * m * v0^2 + 0.5 * K * x0^2;
Em_f = 0.5 * m * v(end)^2 + 0.5 * K * x(end)^2;
deltaEm = Em_f - Em_i;

disp(['Como a variação da Em é diferente de zero (Delta Em = ',num2str(deltaEm),') podemos concluir que não há conservação da Energia'])


%% b)
clear all
close all
clc

K = 1; %N/m
m = 2; %Kg
x0 = 1; %m
v0 = 0; %m/s
w = sqrt(K/m);

h = [2E-1, 1E-1,5E-2, 2E-2,1E-2, 5E-3,2E-3,1E-3]; %s
erro = nan(length(h),1);
x_exact_tf = x0*cos(w*100);

for i = 1:length(h)
    hs = h(i);
    t = 0:hs:100;
    N = length(t);

    x = nan(1,N); x(1) = x0;
    v = nan(1,N); v(1) = v0;

    fv = @(X) -K/m * X;
    fx = @(V) V;

    for k=1:N-1

        r1x = fx(v(k));
        r1v = fv(x(k));

        r2x = fx(v(k)+r1v*hs/2);
        r2v = fv(x(k)+r1x*hs/2);

        r3x = fx(v(k)+(2*r2v-r1v)*hs);
        r3v = fv(x(k)+(2*r2x-r1x)*hs);

        x(k+1) = x(k) +1/6 * (r1x + 4*r2x + r3x) * hs;
        v(k+1) = v(k) +1/6 * (r1v + 4*r2v + r3v) * hs;
    end
    erro(i) = abs(x_exact_tf - x(N));
end

figure(1)
plot(t,x), xlabel('t(s)'), ylabel('x(m)')

figure(2)
plot(log10(h),log10(erro),'ko'), lsline
lsline %traçar uma recta sobre a curva que obteve
aux = polyfit(log10(h),log10(erro),1);  %1 corresponde a ordem do polinomio, como é uma reta a ordem = 1
ordem_metodo = aux(1);
disp(['Ordem do metódo: ',num2str(round(ordem_metodo))])

count = 0;

for k = 2:(N-1) % começamos em 2 pois se começarmos em 1 não temos nenhum anterior para comparar (em matlab não ha posição 0)
    if x(k-1) <= x(k) && x(k) >= x(k+1)
        count = count + 1;
        
        aux = lagr(t(k-1:k+1),x(k-1:k+1));
        t_max(count) = aux(1);
        x_max(count) = aux(2);
    end
end

A = mean(x_max); %ampltiude
disp(['Amplitude média: ',num2str(A),' m'])

p = mean(diff(t_max)) ; %periodo
disp(['Período médio: ',num2str(p),' s'])

%% c)

clear all
close all
clc

h = 0.01;
x0 = 1;v0 = 0;
t0=0;tf=10;



% ode45
reltol = 3E-14;
abstol_1 = 1E-13;
abstol_2 = 1E-13;
options = odeset('RelTol',reltol,'AbsTol',[abstol_1 abstol_2]);
[t_ode45,sol] = ode45(@f,[t0 tf],[x0 v0],options);

figure(1)
plot(t_ode45,sol(:,1),'r--'), xlabel('t'), ylabel('x')
figure(2)
plot(t_ode45,sol(:,2),'r--'), xlabel('t'), ylabel('v')


N = length(sol(:,1));
count = 0;
for k = 2:(N-1) % começamos em 2 pois se começarmos em 1 não temos nenhum anterior para comparar (em matlab não ha posição 0)
    if sol(k-1,1) <= sol(k,1) && sol(k,1) >= sol(k+1,1)
        count = count + 1;
        
        aux = lagr(t_ode45(k-1:k+1),sol(k-1:k+1,1));
        t_max(count) = aux(1);
        x_max(count) = aux(2);
    end
end


disp(['Período médio: ',num2str(t_max),' s'])


function derivadas = f(t,solucao)
K = 1;m = 2;
derivadas = zeros(2,1);

% O vetor solucao tem os valores de x e v para o tempo t em que a função é chamada pela rotina ode45.
derivadas(1) = solucao(2);
% Escreva acima a expressão da derivada de x em função de solucao(1) e de solucao(2).

derivadas(2) = -K/m * solucao(1);
% Escreva a acima expressão da derivada de v em função de solucao(1) e de solucao(2). 
% Se achar que torna o programa mais claro, pode incluir na função as linhas,
% x = solucao(1);
% v = solucao(2);
end