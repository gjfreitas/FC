% F�sica Computacional 2018/2019  Tiago Fernandes
% Problema 7.2a e 7.2b

clear all
close all
clc

%vari�veis shooting
B = 0; %u(0) = 0            %no metodo de shooting B � a solu��o conhecida para u(xinicial)
tol = 10^-12;               %usar 10^-6 ou 10^-12 como tolerancia
nshots = 20;

%vari�vel independente
xmax = 30;
h = 0.01;
x = 0:h:xmax;
V(x>2) = 2.*x(x>2) - 4;
V(x<=2) = 0;
nx = length(x);

% u(r) = r*R(r)         R(r) � a fun��o de onda radial  
psi = zeros(1,nx);

%m�todo de shooting
result = zeros(1, nshots);      %matriz em que v�o ser guardados os resultados u(0) do m�todo de shooting
guess = [0.1 0.2];            %estimativa inicial de E

psi(end) = 0;     %u(N) = 0

psi(end-1) = 1e-15;

%shooting
for ishot = 1:nshots
    
    E = guess(ishot);       %valor de energia usado nesta itera��o
    g = 2*(E - V);        %valor de g(r)
    
    flagmax = 0;        %� 0 at� encontrar o segundo ponto em que E=V(x)   -> indmax
    
    ct1 = 1+ h^2/12 .* g;           %usada no primeiro e terceiro termo do m�todo de numerov
    ct2 = 1- 5*h^2/12 .* g;         %usada no segundo termo do m�todo de numeorv
    
    %m�todo de numerov regressivo
    for i = nx-1:-1:2  

        if abs(E-V(i))<0.05 && flagmax == 0
            indmax = i;
            flagmax = 1;
        end
        psi(i-1) = (ct1(i-1))^-1 * (2*ct2(i) * psi(i) - ct1(i+1)*psi(i+1));
    end
    
    result(ishot) = psi(1);         %result � a o primeiro u -> u(0)
    if ishot >= 2           %a partir da segunda itera��o come�amos a avaliar
                            %os results usando o metodo da secante                         
        m = (result(ishot)-result(ishot-1))/(guess(ishot)-guess(ishot-1));  %declive da secante
        guess(ishot+1) = guess(ishot) + (B - result(ishot))/m;        %calculo da nova guess
        erro = abs(guess(ishot+1)-guess(ishot));        %erro � a diferen�a entre a guess anterior e a nova
        if erro < tol         %quando a guess nova for muito pr�xima da anterior, o m�todo terminou
                              %a ultima guess � a solu��o do problema
            break
        end           
    end
end

%normaliza��o
C = trapz(x, abs(psi).^2);
psinorm = psi ./ sqrt(C);

plot(x, psinorm, 'b')
hold on
plot(x(1:indmax+50),V(1:indmax+50), 'g')

E_obtido = guess(end)

%a probabilidade de estar numa regia�o em que E>= V(x) � dada pelo integral
%abaixo
prob = trapz(x(indmax:end),psinorm(indmax:end).^2)
