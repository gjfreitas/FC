%Teste pratico 1 P2 2010
%% 1.
%% a) P�ndulo grav�tico; Euler Cromer
clc
clear all

%Constantes
m=1; %kg
L=1; %m
g=9.81; %m/s

%Condi��es iniciais
theta0=0.5; %rad
v0=0; %m/s

%Passo
h=0.01;

t=0:h:10;
N=length(t);

theta=zeros(N,1);
v=zeros(N,1);
Ec=zeros(N,1);
Ep=zeros(N,1);

theta(1)=theta0;
v(1)=v0;
Ec(1)=0.5*m*v(1)^2;
Ep(1)=m*g*(L-cos(theta(1)));

for k=1:N-1
    
    v(k+1)=v(k)-(g/L)*sin(theta(k))*h;
    theta(k+1)=theta(k)+v(k+1)*h;
    Ec(k+1)=0.5*m*v(k+1)^2;
    Ep(k+1)=m*g*(L-cos(theta(k+1)));
    
end

%Theta em fun�ao do tempo - sinusoide
plot(t,theta,'r.')
grid on
title('�ngulo em fun�ao do tempo')
xlabel('t (s)')
ylabel('theta (rad)')

%Espa�o de fases - circunferencia
figure
plot(theta,v,'b.')
grid on
title('Espa�o de fases')
xlabel('theta (rad)')
ylabel('velocidade angular (rad/s)')

%Energia (cin�tica, potencial e total) - verificar que a energia total se
%mant�m aproximadamente constante
figure
plot(t,Ec,'r.',t,Ep,'b.',t,Ec+Ep,'y.')
grid on
title('Considera��es energ�ticas')
xlabel('t (s)')
legend('Energia cin�tica','Energia potencial','Energia total')


%% b) Calcular o per�odo
l=1;
for n=1:length(v)-1
    
    if v(n)*v(n+1)<0
        indices(l)=n;
        l=l+1;
    end
    
end

for m=1:length(indices)
    
    
    tempos(m)=interp1(v(indices(m)-2:indices(m)+2),t(indices(m)-2:indices(m)+2),0);
    
end

periodo=2*mean(diff(tempos))

%% c)
clc
clear all

%Constantes
m=1; %kg
L=1; %m
g=9.81; %m/s

%Condi��es iniciais
v0=0; %m/s

%Passo
h=0.01;

t=0:h:10;
N=length(t);

theta=zeros(N,1);
v=zeros(N,1);

periodo_teorico=2*pi*sqrt(L/g);

for n=1:7
    
    theta0=n*0.1;
    theta(1)=theta0;
    v(1)=v0;
    
for k=1:N-1
    
    v(k+1)=v(k)-(g/L)*sin(theta(k))*h;
    theta(k+1)=theta(k)+v(k+1)*h;
    
end

l=1;
for p=1:length(v)-1
    
    if v(p)*v(p+1)<0
        indices(l)=p;
        l=l+1;
    end
    
end

for m=1:length(indices)
    
    
    tempos(m)=interp1(v(indices(m)-2:indices(m)+2),t(indices(m)-2:indices(m)+2),0);
    
end

periodo=2*mean(diff(tempos));

desvio(n)=abs(periodo-periodo_teorico);

end

theta=[0.1 0.2 0.3 0.4 0.5 0.6 0.7];
plot(theta,desvio,'b.')
grid on
xlabel('theta0 (rad)')
ylabel('desvio')
title('Pequenas oscila��es - depend�ncia do desvio com theta0')

%% 2.
%% a)
clc
clear all

%Constantes
m=1; %kg
a=1; %N/m^3

%Condi��es iniciais
x0=1; %m
y0=1.5; %m
v0x=1.2; %m/s
v0y=1; %m/s

%Passo
h=0.01;

t=0:h:50;
N=length(t);

r=zeros(N,2); %vetor posi��o
v=zeros(N,2); %vetor velocidade
Ec=zeros(N,1); %energia cin�tica
Ep=zeros(N,1); %energia potencial

r(1,:)=[x0 y0];
v(1,:)=[v0x v0y];
Ec(1)=0.5*m*(v(1))^2;
Ep(1)=(a/4)*(norm(r(1,:)))^4;

%Usando Euler Cromer
for k=1:N-1
    
    v(k+1,:)=v(k,:)-(a/m)*norm(r(k,:))^2*r(k,:)*h;
    r(k+1,:)=r(k,:)+v(k+1,:)*h;
    Ec(k+1)=0.5*m*(norm(v(k+1,:)))^2;
    Ep(k+1)=(a/4)*(norm(r(k+1,:)))^4;
    

end

%Representa��o gr�fica da trajet�ria
plot(r(:,1),r(:,2),'b.')
grid on
xlabel('Coordenada x (m)')
ylabel('Coordenada y (m)')
title('Trajet�ria')

%Verificar que a energia total (cin�tica+potencial) se mant�m constante
figure
plot(t,Ec,'r.',t,Ep,'b.',t,Ec+Ep,'y.')
grid on
xlabel('tempo (s)')
title('Considera��es energ�ticas')
legend('Energia cin�tica','Energia potencial','Energia total')


%% b)

for n=1:N-1
    
    distancia(n)=norm(r(n,:));
    
end

q=1;
for j=2:length(distancia)-1
    
    if distancia(j)>distancia(j-1) && distancia(j)>distancia(j+1) && distancia(j)>0
        indices1(q)=j;
        q=q+1;

    end
end

for s=1:length(indices1)
    
  t_fino=linspace(t(indices1(s)-2),t(indices1(s)+2),100);
    distancia_max(s)=interp1(t,distancia,t_fino,'spline');
    
    
    
    
end




















