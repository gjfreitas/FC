clear all
close all
clc

lbd = 2.2;

xi = 0;
xf = pi;
h = 0.01;
x = xi:h:xf;
n = length(x);

y = zeros(1,n);
vy = zeros(1,n);

%posicao e velocidade inicial
y(1) = 0;
vy(1) = 1;


%metodo de euler-cromer
fy = @(vy) vy;
fvy = @ (x, y, vy) ((x-lbd)*y+1.6*sin(x)*cos(x)*vy)/(1-0.8*sin(x)^2);
%fv = dv/dt

for i = 1:n-1
    vy(i+1) = vy(i) + fvy( x(i), y(i), vy(i)) * h;
    y(i+1) = y(i) + vy(i+1) * h;            %metodo de euler-cromer usa vx(i+1)
end



%grafico y(x)
plot(x,y, 'y')
set(gca,'Color','k')
xlabel 'x(m)'
ylabel 'y(m)'



