% Teste1_3
%% a)
clc
clear all

%Constantes
g=9.8; %m/s
l=1; %m
b=1; %s^-1

%Condiçoes iniciais
theta0=0.2; %rad
w0=0; %rad/s

%Passo
h=0.01;

t=0:h:30;
N=length(t);

theta=zeros(N,1); %angulo
w=zeros(N,1); %velocidade angular (dtheta/dt)

theta(1)=theta0;
w(1)=w0;

ftheta=@(w) w;
fw=@(g,l,theta,b,w) -(g/l)*sin(theta)-b*w;

for k=1:N-1
    
    r1theta=ftheta(w(k));
    r1w=fw(g,l,theta(k),b,w(k));
    
    r2theta=ftheta(w(k)+r1w*h/2);
    r2w=fw(g,l,theta(k)+r1theta*h/2,b,w(k)+r1w*h/2);
    
    w(k+1)=w(k)+r2w*h;
    theta(k+1)=theta(k)+r2theta*h;
end

plot(t,theta,'r.')
grid on

%% b)
clc
clear all

%Constantes
g=9.8; %m/s
l=1; %m

%Condiçoes iniciais
theta0=0.2; %rad
w0=0; %rad/s

%Passo
h=0.01;

t=0:h:30;
N=length(t);

theta=zeros(N,1); %angulo
w=zeros(N,1); %velocidade angular (dtheta/dt)

theta(1)=theta0;
w(1)=w0;

ftheta=@(w) w;
fw=@(g,l,theta,b,w) -(g/l)*sin(theta)-b*w;

for n=1:10
    
    b=n*0.1;
    
for k=1:N-1
    
    r1theta=ftheta(w(k));
    r1w=fw(g,l,theta(k),b,w(k));
    
    r2theta=ftheta(w(k)+r1w*h/2);
    r2w=fw(g,l,theta(k)+r1theta*h/2,b,w(k)+r1w*h/2);
    
    w(k+1)=w(k)+r2w*h;
    theta(k+1)=theta(k)+r2theta*h;
end

for m=2:length(theta)
    
    if theta(m)>theta(m-1) && theta(m)>theta(m+1) && abs(theta(m))<(theta0/exp(1))
    a=m;
        break
        
    end
    
end
tempo(n)=t(a); 

end

%% c)
clc
clear all

%Constantes
g=9.8; %m/s
l=1; %m
b=1; %s^-1
omega=1;

%Condiçoes iniciais
theta0=0.2; %rad
w0=0; %rad/s

%Passo
h=0.01;

t=0:h:120;
N=length(t);

theta=zeros(N,1); %angulo
w=zeros(N,1); %velocidade angular (dtheta/dt)

theta(1)=theta0;
w(1)=w0;

ftheta=@(w) w;
fw=@(g,l,theta,b,w,omega,t) -(g/l)*sin(theta)-b*w+10*sin(omega*t);

for k=1:N-1
    
    r1theta=ftheta(w(k));
    r1w=fw(g,l,theta(k),b,w(k),omega,t(k));
    
    r2theta=ftheta(w(k)+r1w*h/2);
    r2w=fw(g,l,theta(k)+r1theta*h/2,b,w(k)+r1w*h/2,omega,t(k)+h/2);
    
    w(k+1)=w(k)+r2w*h;
    theta(k+1)=theta(k)+r2theta*h;
end

plot(t,theta,'r.')
grid on













