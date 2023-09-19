clear all
close all
clc

xi = 0;
xf = pi;
h = 0.01;
x = xi:h:xf;
n = length(x);

y = zeros(1,n);
vy = zeros(1,n);

nshots = 100;
B = 0; %y(L) = 0            %no metodo de shooting B é a solução conhecida para y(xfinal)
tol = 10^-12;               %usar 10^-6 ou 10^-12 como tolerancia

result = zeros(1,nshots);        %result contém o valor de y(L) obtido para cada conjunto de guesses

%guess = [0 0.1]        %primeiro valor proprio
guess = [3 3.1];        %segundo valor proprio

for ishot = 1:nshots           %maximo de 100 iterações (o metodo termina muito antes)
    
    y = zeros(1,n);
    vy = zeros(1,n);
    %posicao e velocidade inicial
    y(1) = 0;
    vy(1) = 1;
    
    lbd = guess(ishot);
    fy = @(vy) vy;
    fvy = @ (x, y, vy) ((x-lbd)*y+1.6*sin(x)*cos(x)*vy)/(1-0.8*sin(x)^2);
    
    for i = 1:n-1
        vy(i+1) = vy(i) + fvy( x(i), y(i), vy(i)) * h;
        y(i+1) = y(i) + vy(i+1) * h;            %metodo de euler-cromer usa vx(i+1)
    end

    result(ishot) = y(i+1);         %result é a o ultimo y -> y(L)
    if ishot >= 2           %a partir da segunda iteração começamos a avaliar
                            %os results usando o metodo da secante
                            
        m = (result(ishot)-result(ishot-1))/(guess(ishot)-guess(ishot-1));  %declive da secante
        guess(ishot+1) = guess(ishot) + (B - result(ishot))/m;        %calculo da nova guess
        erro = abs(guess(ishot+1)-guess(ishot));        %erro é a diferença entre a guess anterior e a nova
        if erro < tol         %quando a guess nova for muito próxima da anterior, o método terminou
                              %a ultima guess é a solução do problema
            break
        end
            
    end
end

guess(end)

%grafico y(x)
plot(x,y, 'y')
set(gca,'Color','k')
xlabel 'x(m)'
ylabel 'y(m)'