%Teste1_4
%%1. Decaimento radioativo
%% a)
clc
clear all

%Constantes
lambda1=1;
lambda2=0.2;
lambda3=1;
lambda4=0.4;

%Condições iniciais
c01=1;
c02=0;
c03=0;
c04=0;
c05=0;

c1(1)=c01;
c2(1)=c02;
c3(1)=c03;
c4(1)=c04;
c5(1)=c05;
t(1)=0;

h=0.001;
N=100000;

f1=@(c1,lambda1) -lambda1*c1;
f2=@(c1,lambda1,c2,lambda2) lambda1*c1-lambda2*c2;
f3=@(c2,lambda2,c3,lambda3) lambda2*c2-lambda3*c3;
f4=@(c3,lambda3,c4,lambda4) lambda3*c3-lambda4*c4;
f5=@(c4,lambda4) lambda4*c4;

p=1e-4;

    for k=1:N
        
        t(k+1)=t(k)+h;
        
        r1c1=f1(c1(k),lambda1);
        r1c2=f2(c1(k),lambda1,c2(k),lambda2);
        r1c3=f3(c2(k),lambda2,c3(k),lambda3);
        r1c4=f4(c3(k),lambda3,c4(k),lambda4);
        r1c5=f5(c4(k),lambda4);
        
        r2c1=f1(c1(k)+r1c1*h/2,lambda1);
        r2c2=f2(c1(k)+r1c1*h/2,lambda1,c2(k)+r1c2*h/2,lambda2);
        r2c3=f3(c2(k)+r1c2*h/2,lambda2,c3(k)+r1c3*h/2,lambda3);
        r2c4=f4(c3(k)+r1c3*h/2,lambda3,c4(k)+r1c4*h/2,lambda4);
        r2c5=f5(c4(k)+r1c4*h/2,lambda4);
        
        c1(k+1)=c1(k)+r2c1*h;
        c2(k+1)=c2(k)+r2c2*h;
        c3(k+1)=c3(k)+r2c3*h;
        c4(k+1)=c4(k)+r2c4*h;
        c5(k+1)=c5(k)+r2c5*h;
       
        if c1(k+1)<p && c2(k+1)<p && c3(k+1)<p && c4(k+1)<p
            break
            
        end
        
        
    end
    
 plot(t,c1,'r.',t,c2,'b.',t,c3,'y.',t,c4,'k.',t,c5,'g.')
 grid on
 legend('c1','c2','c3','c4','c5')   
 figure
 
 plot(t,c1+c2+c3+c4+c5,'b.')
 grid on
 
 %% b)
 
 C=0.99*c01;
 l=1;
 for n=1:length(c5)
     
     if c5(n)<C && c5(n+1)>C
         indices1(l)=n;
         l=l+1;
         
         
     end
    
     
 end
 
 
 tempo=interp1(c5(indices1-2:indices1+2),t(indices1-2:indices1+2),C)
 
 %% c)
 s=1;
 q=1;
 z=1;
 
 for n=2:length(c2)-1
     
     if c2(n)>c2(n-1) && c2(n)>c2(n+1)
         
         indices2(s)=n;
         s=s+1;
         
     end
     
     if c3(n)>c3(n-1) && c3(n)>c3(n+1)
     
        indices3(q)=n;
        q=q+1;
        
     end
     
     if c4(n)>c4(n-1) && c4(n)>c4(n+1)
         
         indices4(z)=n;
         z=z+1;
         
     end
     
 end
 
tempo2=t(indices2)
tempo3=t(indices3)
tempo4=t(indices4)

%% d)
clc
clear all

%Constantes
lambda1=1;
lambda2=0.2;
lambda3=1;
lambda4=0.4;

%Condições iniciais
c01=1;
c02=0;
c03=0;
c04=0;
c05=0;

c1(1)=c01;
c2(1)=c02;
c3(1)=c03;
c4(1)=c04;
c5(1)=c05;

N=100000;

f1=@(c1,lambda1) -lambda1*c1;
f2=@(c1,lambda1,c2,lambda2) lambda1*c1-lambda2*c2;
f3=@(c2,lambda2,c3,lambda3) lambda2*c2-lambda3*c3;
f4=@(c3,lambda3,c4,lambda4) lambda3*c3-lambda4*c4;
f5=@(c4,lambda4) lambda4*c4;

p=1e-4;

for n=1:10
    
    h(n)=0.1*n;
    
    for k=1:N
        
        t(1)=0;
        t(k+1,n)=t(k)+h(n);
        
        r1c1=f1(c1(k),lambda1);
%         r1c2=f2(c1(k),lambda1,c2(k),lambda2);
%         r1c3=f3(c2(k),lambda2,c3(k),lambda3);
%         r1c4=f4(c3(k),lambda3,c4(k),lambda4);
%         r1c5=f5(c4(k),lambda4);
        
        r2c1=f1(c1(k)+r1c1*h(n)/2,lambda1);
%         r2c2=f2(c1(k)+r1c1*h(n)/2,lambda1,c2(k)+r1c2*h(n)/2,lambda2);
%         r2c3=f3(c2(k)+r1c2*h(n)/2,lambda2,c3(k)+r1c3*h(n)/2,lambda3);
%         r2c4=f4(c3(k)+r1c3*h(n)/2,lambda3,c4(k)+r1c4*h(n)/2,lambda4);
%         r2c5=f5(c4(k)+r1c4*h(n)/2,lambda4);
        
        c1(k+1,n)=c1(k)+r2c1*h(n);
%         c2(k+1,n)=c2(k)+r2c2*h(n);
%         c3(k+1,n)=c3(k)+r2c3*h(n);
%         c4(k+1,n)=c4(k)+r2c4*h(n);
%         c5(k+1,n)=c5(k)+r2c5*h(n);
       
        if c1(k+1)<p
            break
            
        end
%        indice(1)=0;
%         for m=2:length(c1(:,n))
%             
%            if c1(m,n)==0
%                indice(1)=m;
%                break
%            end
%             
%         end
        
        
    end
    
    c1_analitico=c01*exp(-lambda1*t(end));
    c1_numerico=c1(end);
    erro(n)=abs((c1_numerico)-(c1_analitico));
    
end



plot(h,erro,'r.')
grid on
