%% (1)

%% a)
close all
clear all
clc

epsilon=5;L=0.1;R=10;C=1e-3;
a=1/R/C+R/L;b=2/L/C;c=epsilon/L/C;

t0=0;h=0.000001;tf=10*R*C;
t=t0:h:tf;
N=length(t);
dV=zeros(1,N);V=zeros(1,N);
V(1) = 0;dV(1) = 0;
for i=1:N-1
    dV(i+1)=dV(i)+(-a*dV(i)-b*V(i)+c)*h;
    V(i+1)=V(i)+dV(i)*h;
    if abs(V(i+1)-V(i))<1e-6 && abs(V(i)-epsilon/2)<1e-6
        break
    end 
end
t=t(1:i);V=V(1:i);
plot(t,V,'.')

for j=1:length(V)-1
    if V(j+1)<V(j) && V(j-1)<V(i)
        max=lagr(t(j-1:j+1),V(j-1:j+1));
    end
end
Vmax=max(2)

%% b)
close all
clear all
clc

epsilon=5;L=0.1;R=10;C=1e-3;
a=1/R/C+R/L;b=2/L/C;c=epsilon/L/C;

t0=0;h=0.000001;tf=10*R*C;
t=t0:h:tf;
N=length(t);
dV=zeros(1,N);V=zeros(1,N);
for i=1:N-1
    A=inv([0 1 ; 1 0])*[b*h 1+a*h; 1 -h];
    B=[V(i); dV(i)+h*c];
    aux=linsolve(A,B);
    V(i+1)=aux(1);dV(i+1)=aux(2);
    if abs(V(i+1)-V(i))<1e-6 && abs(V(i)-epsilon/2)<1e-6
        break
    end 
end
t=t(1:i);V=V(1:i);
plot(t,V,'.')
for j=1:length(V)-1
    if V(j+1)<V(j) && V(j-1)<V(i)
        max=lagr(t(j-1:j+1),V(j-1:j+1));
    end
end
Vmax=max(2)

%% c)
clc
close all
clear all

epsilon=5;L=0.1;R=10;C=1e-3;
a=1/R/C+R/L;b=2/L/C;c=epsilon/L/C;

t0=0;h=0.000001;tf=10*R*C;
t=t0:h:tf;
N=length(t);
dV=zeros(1,N);V=zeros(1,N);

fv=@(dV,V) -a*dV-b*V+c;
fx=@(dV) dV;

for i=1:N-1
    r1v=fv(dV(i),V(i));
    r1x=fx(dV(i));  
    r2v=fv(dV(i)+r1x*h/2,V(i)+r1v*h/2);
    r2x=fx(dV(i)+r1v*h/2);
    dV(i+1)=dV(i)+(0*r1v+1*r2v)*h;
    V(i+1)=V(i)+(0*r1x+1*r2x)*h;
    if abs(V(i+1)-V(i))<1e-6 && abs(V(i)-epsilon/2)<1e-6
        break
    end 
end
t=t(1:i);V=V(1:i);dV=dV(1:i);
plot(t,V,'.')


%% (2)

%% a)
close all
clear all
clc

GmS=4*pi^2;x0=1.0167;vy0=8.2;
t0=0;tf=10;h=0.001;
t=t0:h:tf;
N=length(t);
x=zeros(2,N); v=zeros(2,N);theta=zeros(1,N); A=zeros(1,N);
x(1,1)=x0;v(2,1)=vy0;
for i=1:N-1
    r=norm(x(:,i));
    F=-GmS/r^3*x(:,i);
    v(:,i+1)=v(:,i)+h*F;
    x(:,i+1)=x(:,i)+h*v(:,i+1);
    theta(i+1)=mod(atan2(x(2,i),x(1,i)),2*pi);
    A(i)= (((norm(x(:,i+1))+norm(x(:,i)))/2)^2*(theta(i+1)-theta(i)))/2;
    if i~=1 && theta(i+1)<theta(i)
        theta(i+1)=theta(i+1)+2*pi;
        break
    end
end

T=interp1(theta(i-1:i+1),t(i-1:i+1),2*pi)
x=x(:,1:i);v=v(:,1:i);t=t(1:i);theta=theta(1:i);A=A(1:i);

figure(1)
plot(0,0,'oy','MarkerSize',10,'Linewidth',20);
hold on
plot(x(1,:),x(2,:),'.')
legend('Sol','órbita da Terra')
axis equal
% set(gca,'PlotBoxAspectRatio',[1 1 1])


