%% 1o Teste Prático de Física Computacional — 2010/2011 17 de Março de 2011
clear all
close all
clc

%a=G*ms
a=4*pi^2;

x0 = 49.305; % Unidades Astronomicas (AU)
y0 = 0; % AU
v0x = 0; % AU/ano
v0y = 0.775566; % AU/ano

% Metódo de Euler-Cromer

h=0.0001;
t=0:h:500;
N=length(t);

r=nan(N,2);
v=nan(N,2);
Ec=nan(N,1);
Ep=nan(N,1);
Er=nan(N,1);

r(1,:)=[x0 y0];
v(1,:)=[v0x v0y];
Ec(1)=0.5*norm(v(1,:))^2;
Ep(1)=-a/norm(r(1,:));

% a = G * ms = 4 * pi^2
% Ec/Ep = (1/2 * mc * v^2)/(-a*m/r) = (-v^2 * r)/ 2a
Er(1)=-(norm(v(1,:))^2*norm(r(1,:)))/(2*a);


for k=1:N-1
    
    v(k+1,:)=v(k,:)-(a/(norm(r(k,:)))^3)*r(k,:)*h;
    r(k+1,:)=r(k,:)+v(k+1,:)*h;
    Ec(k+1)=0.5*(norm(v(k+1,:)))^2;
    Ep(k+1)=-a/norm(r(k+1,:));
    Er(k+1)=-((norm(v(k+1,:))^2)*norm(r(k+1,:)))/(2*a);
    
end

plot(r(:,1),r(:,2),'r.'), xlabel('x (AU)'), ylabel('v (AU/ano)')
grid on

% d)
figure
plot(t,Ec,'r.',t,Ep,'b.',t,Er,'m.'), xlabel('t(ano)'), ylabel('E(J)')
grid on
legend('Ec','Ep','Ec/Ep')
Em = Ec+Ep;
figure
plot(t,Em), xlabel('t (ano)'), ylabel('Em (J)')
grid on

% b) Encontrar o perido

x=r(:,1);
y=r(:,2);

l=1;
for n=1:length(y)-1
    if y(n)*y(n+1)<0
        indices(l)=n;
        l=l+1;
        
    end
end

for m=1:length(indices)
    
    tempos(m)=interp1(y(indices(m)-2:indices(m)+2),t(indices(m)-2:indices(m)+2),0);
    
end
periodo=2*mean(diff(tempos));

disp(['periodo: ',num2str(periodo),' anos'])

% c) Determinar os semi-eixos da órbita

vx=v(:,1);
vy=v(:,2);

%Encontrar o eixo vertical
p=1;
for j=1:length(vy)-1
    
    if vy(j)*vy(j+1)<0 && y(j)>0
        indices2(p)=j;
        p=p+1;
        
    end

end

%Encontrar a distancia minima ao Sol para calcular o eixo horizontal
q=1;
for i=1:length(y)-1
    
    if y(i)*y(i+1)<0 && x(i)<0
        indices3(q)=i;
        q=q+1;
    end
    
end

semieixo_horizontal=(abs(mean(x(indices3)))+x0)/2;
semieixo_vertical=mean(y(indices2));
disp(['Semi Eixo Horizontal : ',num2str(semieixo_horizontal),' AU'])
disp(['Semi Eixo Vertical : ',num2str(semieixo_vertical),' AU'])

% e) Velocidade média e perímetro da órbita

for s=1:N
    
    vmod(s)=norm(v(s,:));
    
end

v_media=(mean(vmod)*150E6)/(365.25*24*3600);
disp(['Velocidade média : ',num2str(v_media),' m/s'])

theta=mod(atan2(y,x),2*pi);

t=0:h:periodo;
P=0;
for k=1:length(t)-2
    
    P=P+norm(r(k,:))*(theta(k+1)-theta(k));
    
end
disp(['Perimetro = ',num2str(P),' UA'])

% f)

t=0:h:periodo/4;
theta=mod(atan2(y,x),2*pi);
A1=0;
for n=1:length(t)-1
    
    A1=A1+(norm(r(n,:))^2)*(theta(n+1)-theta(n))/2;
    
    
end

t=periodo/4:h:periodo/2;
N=length(t);
A2=0;

for m=1:length(t)-1
    
    A2=A2+(norm(r(n,:))^2)*(theta(n+1)-theta(n))/2;
    
end

disp(['A1 = ',num2str(A1),' UA ^ 2'])
disp(['A2 = ',num2str(A2),' UA ^ 2'])
disp('A1 = A2, Logo a Lei de Kepler confirma-se')
