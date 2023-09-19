%Teste1_5
%% 1. Osciladores acoplados
%% a)

%Constantes
m1=1;
m2=2;
k1=1;
k2=0.8;
r1=2;
r2=2;

%Condiçoes iniciais
x01=1.5;
x02=4.2;
v01=0;
v02=0;

%Passo
h=0.001;

t=0:h:80;
N=length(t);

x1=zeros(N,1);
v1=zeros(N,1);
x1(1)=x01;
v1(1)=v01;

x2=zeros(N,1);
v2=zeros(N,1);
x2(1)=x02;
v2(1)=v02;

fx1=@(v1) v1;
fv1=@(m1,x1,r1,k1,x2,r2,k2) (1/m1)*(-k1*(x1-r1)+k2*(x2-x1-r2));

fx2=@(v2) v2;
fv2=@(m2,k2,x1,x2,r2) -(k2/m2)*(x2-x1-r2);


for k=1:N-1
    
    r1x1=fx1(v1(k));
    r1x2=fx2(v2(k));
    r1v1=fv1(m1,x1(k),r1,k1,x2(k),r2,k2);
    r1v2=fv2(m2,k2,x1(k),x2(k),r2);
    
    r2x1=fx1(v1(k)+r1v1*h/2);
    r2x2=fx2(v2(k)+r1v2*h/2);
    r2v1=fv1(m1,x1(k)+r1x1*h/2,r1,k1,x2(k)+r1x2*h/2,r2,k2);
    r2v2=fv2(m2,k2,x1(k)+r1x1*h/2,x2(k)+r1x2*h/2,r2);
    
    x1(k+1)=x1(k)+r2x1*h;
    x2(k+1)=x1(k)+r2x2*h;
    v1(k+1)=v1(k)+r2v1*h;
    v2(k+1)=v2(k)+r2v2*h;
    
end

plot(t,x1,'b.',t,x2,'r.')
grid on
title('Posição dos dois corpos ao longo do tempo')
xlabel('tempo(s)')
ylabel('x(cm)')
legend('Corpo 1','Corpo 2')

%% b)

l=1;
for m=1:length(x2)-1
    
   if x2(m)<4 && x2(m+1)>4 && v2(m)>0
    indices1(l)=m;
    l=l+1;
    
   end
   
end




