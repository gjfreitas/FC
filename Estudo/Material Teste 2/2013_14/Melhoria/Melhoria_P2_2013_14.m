%% Melhoria 2o Teste Prático Física Computacional — 2013/2014
%% 1 a)
clear all
close all
clc

% Constantes
m = 1;
R = 8E-2;
p = 1.22;
C = 0.47;
g = 9.8;
const = (pi * R^2)/(2*m) * C;

% Condições iniciais
x0 = 0;
y0 = 0;
ang = 30;
v0 = 100;

v0x = v0*cosd(ang);
v0y = v0*sind(ang);

 t0 = 0; %s
 tf = 1000; %s
 h = 0.01; %s
 t = t0:h:tf; %s
 
 N = length(t);
 x = nan(1,N); x(1) = x0;
 y = nan(1,N); y(1) = y0;
 vx = nan(1,N); vx(1) = v0x;
 vy = nan(1,N);vy(1) = v0y;
 
for k = 1:N-1
     if y(k) < 0
         break %penultimo valor de z positivo e ultimo negativo, isto faz-se para se poder usar interp1
     end

    vx(k+1) = vx(k) + (-const *sqrt((vy(k)^2 + vx(k)^2)) * vx(k))*h;
    x(k+1) = x(k) + vx(k+1)*h;

    vy(k+1) = vy(k) + (-const * sqrt((vy(k)^2 + vx(k)^2)) * vy(k) - g)*h;
    y(k+1) = y(k)+ vy(k+1)*h;

end
 
plot(x,y,'m-'), xlabel('x'), ylabel('y'), title('Trajetória')

t_impacto = interp1([y(k-1),y(k)],[t(k-1),t(k)],0);
x_impacto = interp1([t(k-1),t(k)],[x(k-1),x(k)],t_impacto);
disp(['tempo de voo: ',num2str(t_impacto),' s'])
disp(['alcance: ',num2str(x_impacto),' m'])

%% 1 b)
clear all
close all
clc

% Constantes
m = 1;
R = 8E-2;
p = 1.22;
C = 0.47;
g = 9.8;

const = (pi * R^2)/(2*m) * C;

% Condições iniciais
x0 = 0;
y0 = 0;
ang = 30;
v0 = 100;

t0 = 0; %s
tf = 1000; %s
h = 0.01; %s
t = t0:h:tf; %s

N = length(t);
x = nan(1,N); x(1) = x0;
y = nan(1,N); y(1) = y0;

% shooting
B = 240.0;
guess = [50.0 54.0];
tol = 1E-5;
nshots = 1000;

for j = 1:nshots
    clear ang
    clear v0y
    clear v0x
    
    
    ang = guess(j);
    v0x = v0*cosd(ang);
    v0y = v0*sind(ang);
    vx = nan(1,N); vx(1) = v0x;
    vy = nan(1,N);vy(1) = v0y;

    for k = 1:N-1
        if y(k) < 0
            break %penultimo valor de z positivo e ultimo negativo, isto faz-se para se poder usar interp1
        end
        vx(k+1) = vx(k) + (-const*sqrt((vy(k)^2 + vx(k)^2)) * vx(k))*h;
        x(k+1) = x(k) + vx(k+1)*h;

        vy(k+1) = vy(k) + (-const * sqrt((vy(k)^2 + vx(k)^2)) * vy(k) - g)*h;
        y(k+1) = y(k)+ vy(k+1)*h;
    end

    result(j) = x(k);
    diff = B - result(j);
    
    if j>= 2
        m = (result(j)-result(j-1))/(guess(j)-guess(j-1));
        guess(j+1) = guess(j)+(diff)/m;
        if abs(guess(j)-guess(j-1)) < tol
            break
        end
    end
    
end

disp(['alcance: ',num2str(result(j)),' m'])
disp(['ang: ',num2str(guess(j)),' graus'])
plot(x,y,'m-'), xlabel('x'), ylabel('y'), title('Trajetória')

%% 2
clear all
close all
clc

% No enuncidado n falam de valor de lamba então vou considerar = 1
lambda = 1;
h = 0.01;
x = -5:h:5;
N = length(x);
const = h^2 * lambda;

% Construção da matriz A
for i = 1:N-2
    A1(i,1) = -2-h^2*(x(i)^2);
end
A1 = diag(A1); % transforma a matriz A1 numa matriz diagonal
A2 = diag(ones([1 N-3]),1); % sobe um posição relativamente a diagonal
A3 = diag(ones([1 N-3]),-1);% desce uma posição relativamente a diagonal
A = A1+A2+A3;

% phi = y = > Ay = const y (ver ultimo slide - Aula 6 e trabalho 5.2)

% calcular os vetores e valores próprios de A
[vec,val]=eigs(A,3,'sm');

% Para calcular os valores próprios de phi tem que se dividir por -h^2 * lambda
v1 = val(1,1)/const;
v2 = val(2,2)/const;
v3 = val(3,3)/const;

disp(['Valor próprio 1 ~ ',num2str(round(v1))]) 
disp(['Valor próprio 2 ~ ',num2str(round(v2))])
disp(['Valor próprio 3 ~ ',num2str(round(v3))])



plot(x(2:(end-1)),vec(:,1),'r-',x(2:(end-1)),vec(:,2),'m-',x(2:(end-1)),vec(:,3),'b-')
xlabel('x'),ylabel('vetor próprio')
legend('Vetor próprio 1', 'Vetor próprio 2', 'Vetor próprio 3')



