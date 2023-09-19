%Teste pratico 1 2012 P3
%% a)
clc
clear all

m=1.5;
K=2;
a=-0.1;
x0=1.9;
v0=0.1;

h=0.01;
t=0:h:30;
N=length(t);

x=zeros(N,1);
v=zeros(N,1);
E=zeros(N,1);
x(1)=x0;
v(1)=v0;
E(1)=0.5*m*v(1)^2+(0.5*K*x(1)^2)*(1+a*x(1));

for k=1:N-1
    
    v(k+1) = v(k)-(K/m)*x(k)*(1+(3/2)*a*x(k))*h;
    x(k+1) = x(k)+v(k+1)*h;
    E(k+1) = 0.5*m*v(k+1)^2+(0.5*K*x(k+1)^2)*(1+a*x(k+1));

    
end

plot(x,v,'r.')
grid on
figure
plot(t,x,'b.')
grid on
figure
plot(t,E,'y.')
grid on

m=1;
for n=1:length(x)-1
    if x(n+1)*x(n)<0
       indices1(m)=n; 
       m=m+1;
    end
    
end

u=1;
for j=1:length(indices1)
    
    p(u)=interp1(x(indices1(j)-2:indices1(j)+2),t(indices1(j)-2:indices1(j)+2),0);
    u=u+1;
  
end

periodo=2*mean(diff(p))

i=1;
q=1;
for q=2:length(x)-1
    if x(q)>x(q+1) && x(q)>x(q-1) && x(q)>0
        indices2(i)=q;
        i=i+1;
        
    end
    
end
    

s=1;
for r=1:length(indices2)
    
    tfino=t(indices2(r)-2):0.001:t(indices2(r)+2);
    
    amp=interp1(t,x,tfino);
    A(s)=mean(amp);
    s=s+1;

end

amplitude=mean(A)

%% b)
clc
clear all

m=1.5;
K=2;
a=-0.1;
x0=0;
v0=0.1;

h=0.01;
t=0:h:30;
N=length(t);

x=zeros(N,1);
v=zeros(N,1);
E=zeros(N,1);
x(1)=x0;
v(1)=v0;
E(1)=0.5*m*v(1)^2+(0.5*K*x(1)^2)*(1+a*x(1));










