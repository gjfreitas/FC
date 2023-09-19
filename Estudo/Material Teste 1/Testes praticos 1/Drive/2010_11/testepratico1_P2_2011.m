% Teste pratico 1 P2 2011
%% a) Usando Euler Cromer

clc
clear all

%a=G*ms
a=4*pi^2;
x0=49.305;
y0=0;
v0x=0;
v0y=0.775566;

h=0.0001;
t=0:h:500;
N=length(t);

r=zeros(N,2);
v=zeros(N,2);
Ec=zeros(N,1);
Ep=zeros(N,1);

r(1,:)=[x0 y0];
v(1,:)=[v0x v0y];
Ec(1)=0.5*norm(v(1,:))^2;
Ep(1)=-a/norm(r(1,:));

for k=1:N-1
    
    v(k+1,:)=v(k,:)-(a/(norm(r(k,:)))^3)*r(k,:)*h;
    r(k+1,:)=r(k,:)+v(k+1,:)*h;
    Ec(k+1)=0.5*(norm(v(k+1,:)))^2;
    Ep(k+1)=-a/norm(r(k+1,:));
    
end

plot(r(:,1),r(:,2),'r.')
grid on
figure
plot(t,Ec,'r.',t,Ep,'b.',t,Ec+Ep,'y.')
grid on


%% b) Calcular o período
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

periodo=2*mean(diff(tempos))

%% c) Determinar os semi-eixos da órbita

vx=v(:,1);
vy=v(:,2);

%Encontrar o eixo vertical
p=1;
for j=1:length(vy)-1;
    
    if vy(j)*vy(j+1)<0 && y(j)>0
        indices1(p)=j;
        p=p+1;
        
    end

end

%Encontrar a distancia minima ao Sol para calcular o eixo horizontal
q=1;
for i=1:length(y)-1;
    
    if y(i)*y(i+1)<0 && x(i)<0
        indices2(q)=i;
        q=q+1;
    end
    
end

semieixo_horizontal=(abs(mean(x(indices2)))+x0)/2
semieixo_vertical=mean(y(indices1))

%% e) Velocidade média e perímetro da órbita

for s=1:N
    
    vmod(s)=norm(v(s,:));
    
end

v_media=(mean(vmod)*150e6)/(365.25*24*3600)

theta=mod(atan2(y,x),2*pi);

t=0:h:periodo;
P=0;
for k=1:length(t)-2
    
    P=P+norm(r(k,:))*(theta(k+1)-theta(k));
    
end

%% f)

t=0:h:periodo/4;
N=length(t);

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