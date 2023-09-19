%Teste prático 2012 P1
%1
%% a)
clc
clear all

%a=G*ms
x0=-0.1246177;
y0=0;
v0x=0;
v0y=-24.9099;
r0=[x0 y0];
v0=[v0x v0y];
a=4*pi^2;
h=0.0001;
t=0:h:10;
N=length(t);
r=zeros(N,2);
v=zeros(N,2);
r(1,:)=r0;
v(1,:)=v0;

for k=1:N-1
    
    v(k+1,:)=v(k,:)-(a/(norm(r(k,:)))^3)*r(k,:)*h;
    r(k+1,:)=r(k,:)+v(k+1,:)*h;
    
end

plot(r(:,1),r(:,2),'r.')
grid on
axis([-6 6 -6 6])
set(gca,'PlotBoxAspectRatio',[1 1 1])

%b)

x=r(:,1);
y=r(:,2);

l=1;
for j=1:length(y)-1
    
    if y(j)*y(j+1)<0 && x(j)>0
    indices(l)=j;
    l=l+1;
    end
    
end

for q=1:length(indices)
    
    p(q)=interp1(y(indices(q)-2:indices(q)+2),t(indices(q)-2:indices(q)+2),0,'spline');
      
end

periodo=mean(diff(p))


for m=1:length(indices)
    
    d(m)=interp1(y(indices(m)-2:indices(m)+2),x(indices(m)-2:indices(m)+2),0,'spline');
   
   
end

distancia_max=mean(d)

%c)
T=0:h:periodo;
theta=mod(atan2(r(:,2),r(:,1)),2*pi);
A=0;
dtheta=diff(theta);

for u=1:length(T)-1
    
    A=A+(((norm(r(u,:))^2)*(dtheta(u)))/2); 
end

Area=A


%d)

x0=-0.1246177;
y0=0;
v0x=0;
v0y=-24.9099;
r0=[x0 y0];
v0=[v0x v0y];
%a=G*ms
a=4*pi^2;
h=0.0001;
t=0:h:periodo/2;
N=length(t);
r=zeros(N,2);
v=zeros(N,2);
r(1,:)=r0;
v(1,:)=v0;
Er=zeros(N,1);
Er(1)=-((v(1))^2*norm(r(1,:)))/(2*a);

for k=1:N-1
    
    v(k+1,:)=v(k,:)-(a/(norm(r(k,:)))^3)*r(k,:)*h;
    r(k+1,:)=r(k,:)+v(k+1,:)*h;
    Er(k+1)=-((v(k+1))^2*norm(r(k+1,:)))/(2*a);
end

figure
plot(t,Er,'b.')
grid on

%% 2

clc
clear all

x0=-0.9141335;
y0=0;
v0x=0;
v0y=-9.2822965;
r0=[x0 y0];
v0=[v0x v0y];
%a=G*ms

a=4*pi^2;

N=4000000;
r=zeros(N,2);
v=zeros(N,2);
r(1,:)=r0;
v(1,:)=v0;


t(1)=0;
for k=1:N-1
    
    if norm(r(k,:))<=70
    h=0.00005;
    
    else 
    h=0.0025;
    end
    
    t(k+1)=t(k)+h;
    
    v(k+1,:)=v(k,:)-(a/(norm(r(k,:)))^3)*r(k,:)*h;
    r(k+1,:)=r(k,:)+v(k+1,:)*h;
        
    end


plot(r(:,1),r(:,2),'b.')
grid on
axis([-100 400 -400 400])



