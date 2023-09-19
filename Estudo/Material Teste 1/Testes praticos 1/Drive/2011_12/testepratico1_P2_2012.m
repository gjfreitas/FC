%Teste pratico1 P2 2012
%1
%% a)
clc
clear all

g=9.81;
m=0.145;
R=36.6*10^-3;
A=pi*R^2;
x0=0;
y0=0.9;
v0=45;
v0x=45*cosd(30);
v0y=45*sind(30);
ro=1.225;

h=0.01;
t=0:h:10;
N=length(t);

v=zeros(N,2);
r=zeros(N,2);
v(1,:)=[v0x v0y];
r(1,:)=[x0 y0];
P=[0 -m*g];
e=0;
for k=1:N-1
    
    Cd=0.2194+0.3262/(1+exp((norm(v(k,:))-35)/5));
    Fd=(-0.5*Cd*ro*A*norm(v(k,:)))*v(k,:);
    
    v(k+1,:)=v(k,:)+(1/m)*(P+Fd)*h;
    r(k+1,:)=r(k,:)+v(k,:)*h;
    
    e=e+sqrt((r(k,1)-r(k+1,1))^2+(r(k,2)-r(k+1,2))^2);
end

plot(r(:,1),r(:,2),'r.')
grid on

x=r(:,1);
y=r(:,2);

alcance=interp1(y,x,0)
tempo_voo=interp1(y,t,0)

for n=1:length(v)
    
    vescalar(n)=norm(v(n,:));
    
end

vmedia=sum(vescalar)/length(v)

%% b)
clc
clear all

g=9.81;
m=0.145;
R=36.6*10^-3;
A=pi*R^2;
x0=0;
y0=0.9;
v0=45;
v0x=45*cosd(30);
v0y=45*sind(30);
ro=1.225;
f=25/3;

h=0.01;
t=0:h:10;
N=length(t);

v=zeros(N,3);
r=zeros(N,3);
v(1,:)=[v0x v0y 0];
r(1,:)=[x0 y0 0];
P=[0 -m*g 0];
w=[0 0 2*pi*f];

for k=1:N-1
    
    Cd=0.2194+0.3262/(1+exp((norm(v(k,:))-35)/5));
    Fd=(-0.5*Cd*ro*A*norm(v(k,:)))*v(k,:);
    S=(R*norm(w))/norm(v(k,:));
    
    if S<=0.1
        Cl=1.5*S;
    else
        Cl=0.09+0.6*S;
    end
    Fl=0.5*Cl*ro*A*norm(v(k,:))^2*(cross((w/norm(w)),v(k,:)/norm(v(k,:))));
    
    v(k+1,:)=v(k,:)+(1/m)*(P+Fd+Fl)*h;
    r(k+1,:)=r(k,:)+v(k,:)*h;
    
end

x=r(:,1);
y=r(:,2);
plot(x,y,'r.');
grid on

alcance=interp1(y,x,0)

%% c)
clc
clear all

g=9.81;
m=0.145;
R=36.6*10^-3;
A=pi*R^2;
x0=0;
y0=0.9;
ro=1.225;
f=25/3;

h=0.01;
t=0:h:10;
N=length(t);

v=zeros(N,3);
r=zeros(N,3);

r(1,:)=[x0 y0 0];
P=[0 -m*g 0];
w=[0 0 2*pi*f];

for n=1:89
    theta(n)=n;
    v0x=45*cosd(n);
    v0y=45*sind(n);
    v(1,:)=[v0x v0y 0];
    
for k=1:N-1
    
    Cd=0.2194+0.3262/(1+exp((norm(v(k,:))-35)/5));
    Fd=(-0.5*Cd*ro*A*norm(v(k,:)))*v(k,:);
    S=(R*norm(w))/norm(v(k,:));
    
    if S<=0.1
        Cl=1.5*S;
    else
        Cl=0.09+0.6*S;
    end
    Fl=0.5*Cl*ro*A*norm(v(k,:))^2*(cross((w/norm(w)),v(k,:)/norm(v(k,:))));
    
    v(k+1,:)=v(k,:)+(1/m)*(P+Fd+Fl)*h;
    r(k+1,:)=r(k,:)+v(k,:)*h;
    
end

x=r(:,1);
y=r(:,2);

alcance(n)=interp1(y,x,0);

end

plot(theta,alcance,'b*')
grid on
l=1;
for j=2:length(alcance)-1
    
    if alcance(j)>alcance(j-1)&&alcance(j)>alcance(j+1)
    indices(l)=j;
    l=l+1;
    end
end

%interp1(alcance(indices-2:indices+2),theta(indices-2;indices+2),

%d)
%% Calculando analiticamente a derivada
clc
clear all

g=9.81;
m=0.145;
R=36.6*10^-3;
A=pi*R^2;
x0=0;
y0=0.9;
v0=45;
v0x=45*cosd(30);
v0y=45*sind(30);
ro=1.225;
f=25/3;

h=0.01;
t=0:h:10;
N=length(t);

v=zeros(N,3);
r=zeros(N,3);
v(1,:)=[v0x v0y 0];
r(1,:)=[x0 y0 0];
P=[0 -m*g 0];
w=[0 0 2*pi*f];

for k=1:N-1
    
    Cd=0.2194+0.3262/(1+exp((norm(v(k,:))-35)/5));
    dCd=(-(0.3263/5)*exp(norm(v(k,:))/5-7))/(1+exp(norm(v(k,:))/5-7))^2;
    Fd=(-0.5*Cd*ro*A*norm(v(k,:)))*v(k,:);
    S=(R*norm(w))/norm(v(k,:));
    
    Cl=2*S*Cd*(1+(norm(v(k,:))/2*Cd))*dCd;
    Fl=0.5*Cl*ro*A*norm(v(k,:))^2*(cross((w/norm(w)),v(k,:)/norm(v(k,:))));
    
    v(k+1,:)=v(k,:)+(1/m)*(P+Fd+Fl)*h;
    r(k+1,:)=r(k,:)+v(k,:)*h;
    
end

x=r(:,1);
y=r(:,2);
plot(x,y,'r.');
grid on

alcance=interp1(y,x,0)

%% Calculando numericamente a derivada
clc
clear all

g=9.81;
m=0.145;
R=36.6*10^-3;
A=pi*R^2;
x0=0;
y0=0.9;
v0=45;
v0x=45*cosd(30);
v0y=45*sind(30);
ro=1.225;
f=25/3;

h=0.01;
t=0:h:10;
N=length(t);

v=zeros(N,3);
r=zeros(N,3);
v(1,:)=[v0x v0y 0];
r(1,:)=[x0 y0 0];
P=[0 -m*g 0];
w=[0 0 2*pi*f];

for k=2:N-1
    Cd=0.2194+0.3262/(1+exp((norm(v(k,:))-35)/5));
    dCd=(0.2194+0.3262/(1+exp((norm(v(k,:))-35)/5))-0.2194+0.3262/(1+exp((norm(v(k-1,:))-35)/5)))/(norm(v(k,:))-norm(v(k-1,:)));
    Fd=(-0.5*Cd*ro*A*norm(v(k,:)))*v(k,:);
    S=(R*norm(w))/norm(v(k,:));
    
     Cl=2*S*Cd*(1+(norm(v(k,:))/2*Cd))*dCd;
    Fl=0.5*Cl*ro*A*norm(v(k,:))^2*(cross((w/norm(w)),v(k,:)/norm(v(k,:))));
    
    v(k+1,:)=v(k,:)+(1/m)*(P+Fd+Fl)*h;
    r(k+1,:)=r(k,:)+v(k,:)*h;
    
end

x=r((2:end),1);
y=r((2:end),2);
plot(x,y,'r.');
grid on

alcance=interp1(y,x,0)

