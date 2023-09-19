%% 1o Teste Prático Física Computacional — 2013/2014 13 março de 2014
% Turma P4
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

h=0.001;
tf = 1000; % anos
t = 0:h:tf;
N=length(t);

c1 = nan(1,N); c1(1)=c01;
c2 = nan(1,N); c2(1)=c02;
c3 = nan(1,N); c3(1)=c03;
c4 = nan(1,N); c4(1)=c04;
c5 = nan(1,N); c5(1)=c05;



f1=@(c1,lambda1) -lambda1*c1;
f2=@(c1,lambda1,c2,lambda2) lambda1*c1-lambda2*c2;
f3=@(c2,lambda2,c3,lambda3) lambda2*c2-lambda3*c3;
f4=@(c3,lambda3,c4,lambda4) lambda3*c3-lambda4*c4;
f5=@(c4,lambda4) lambda4*c4;

p=1e-4;

for k=1:N

    

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

    if c1(k+1) < p && c2(k+1) < p && c3(k+1) < p && c4(k+1) < p
        break

    end


end
 
 figure(1)
 plot(t,c1,'r.',t,c2,'b.',t,c3,'y.',t,c4,'k.',t,c5,'g.'), xlabel('t(anos)')
 grid on
 legend('c1','c2','c3','c4','c5')  
 
 figure(2)
 plot(t,c1+c2+c3+c4+c5,'b.')
 grid on
 
 % b)

 C=0.99*c01;
 lcount=1;
 for n=1:length(c5)
     
     if c5(n)<C && c5(n+1)>C
         indices1(lcount)=n;
         lcount=lcount+1;
     end
 end
 

 tempo=interp1(c5(indices1-2:indices1+2),t(indices1-2:indices1+2),C);
 disp(['Tempo para o qual a população de c5 é 99% da poupulação c1 : ', num2str(tempo), ' anos'])
 
 %c)
 s=1;
 q=1;
 z=1;
 K = length(c2);
 
 for n=2:K
     
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
 

tempo2 = max(t(indices2));
disp(['c2 tem maximo a : ', num2str(tempo2), ' anos'])
tempo3 = max(t(indices3));
disp(['c3 tem maximo a : ', num2str(tempo3), ' anos'])
tempo4 = max(t(indices4));
disp(['c4 tem maximo a : ', num2str(tempo4), ' anos'])

%% d)
clear all
close all
clc

tf = 100;
lambda1 = 1;
c01 = 1;
c1(1) = c01;

N = 1000;
f1 = @(c1,lambda1) -lambda1*c1;

for n=1:10
    
    h(n) = 0.1*n;
    
    for k=1:N
        
        t(1) = 0;
        t(k+1,n) = t(k)+h(n);

        r1c1 = f1(c1(k),lambda1);
        r2c1 = f1(c1(k)+r1c1*h(n)/2,lambda1);
        c1(k+1,n) = c1(k)+r2c1*h(n);

        if c1(k+1) < 1e-4
            break     
        end   
    end
    
    c1_analitico=c01*exp(-lambda1*t(end));
    c1_numerico=c1(end);
    erro(n)=abs((c1_numerico)-(c1_analitico));
    
end

plot(h,erro,'r.','MarkerSize',10)
grid on


% lsline %traçar uma recta sobre a curva que obteve
% aux = polyfit(log10(h),log10(erro),2);  %2 corresponde a ordem do polinomio, como é uma reta a ordem = 2
% declive = aux(1);
% disp([' declive: ',num2str(declive)])
