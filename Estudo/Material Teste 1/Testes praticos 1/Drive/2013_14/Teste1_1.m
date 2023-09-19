%Teste1_1
%% a)
clc
clear all

%Constantes
g=9.8; %m/s^2
r=0.03; %m

%Condições iniciais
x0=0.025; %m
V0=(pi/3)*(15.625); %m^3
A0=15.625/x0; %m^2

R=sqrt(A0/pi);

theta=atan(R/x0) %rad

%% b)
clc
clear all

%Constantes
g=9.8; %m/s^2
r=0.03; %m

%Condições iniciais
x0=0.025; %m
V0=(pi/3)*(15.625); %m^3
A0=15.625/x0; %m^2

R=sqrt(A0/pi);

theta=atan(R/x0); %rad

%Passo 
h=2;

tf=8*60; %seg

t=0:h:tf;
N=length(t);

x=zeros(N,1);
x(1)=x0;

%A(x)=pi*(tan(theta))^2*x^2

f=@(r,g,x,theta,pi) -0.6*pi*r^2*(sqrt(2*g*x)/(pi*(tan(theta))^2*x^2));

for k=1:N-1
    
    r1=f(r,g,x(k),theta,pi);
    r2=f(r,g,x(k)+r1*h/2,theta,pi);
    r3=f(r,g,x(k)+r2*h/2,theta,pi);
    r4=f(r,g,x(k)+r3*h/2,theta,pi);
    
    x(k+1)=x(k)+(1/6)*(r1+2*r2+2*r3+r4)*h;
    
    if x(k+1)<0
        
        break
        
    end
end

plot(t,x,'r.')
grid on

%% c)
clc
clear all

%Constantes
g=9.8; %m/s^2
r=0.03; %m

%Condições iniciais
x0=0.025; %m
V0=(pi/3)*(15.625); %m^3
A0=15.625/x0; %m^2

R=sqrt(A0/pi);

theta=atan(R/x0); %rad

%Passo 
h=[2, 2.5, 3, 3.5, 4];

tf=8*60; %seg

t=0:h:tf;
N=length(t);

x=zeros(N,1);
x(1)=x0;

%A(x)=pi*(tan(theta))^2*x^2

f=@(r,g,x,theta,pi) -0.6*pi*r^2*(sqrt(2*g*x)/(pi*(tan(theta))^2*x^2));

for n=1:5
    
    h(n)=n;
    t=0:h(n):tf;
    N=length(t);

    x=zeros(N,1);
    x(1)=x0;
    
for k=1:N-1
    
    r1=f(r,g,x(k),theta,pi);
    r2=f(r,g,x(k)+r1*h(n)/2,theta,pi);
    r3=f(r,g,x(k)+r2*h(n)/2,theta,pi);
    r4=f(r,g,x(k)+r3*h(n)/2,theta,pi);
    
    x(k+1,n)=x(k)+(1/6)*(r1+2*r2+2*r3+r4)*h(n);
    
    if x(k+1)<0
        
        break
        
    end
    
    estimativa_num(n)=x(end);
end

end



