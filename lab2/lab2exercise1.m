clc;
clear;

#y' = −100y + 100t + 101
#y(0) = y0.

function [] = lab2exercise1(y0, h, N)

    yEuler = zeros(N+1,1);
    ybEuler = zeros(N+1,1);
    yEuler(1,1) = y0;
    ybEuler(1,1) = y0;

    for n = 1:N
        t = (n-1)*h;
        yEuler(n+1,1) = yEuler(n,1) + h*(-100*yEuler(n,1) + 100*t +101);
        ybEuler(n+1,1) = (ybEuler(n,1) + h*(100*(t+h) +101))/(1 + 100*h);
    end

    figure;
    plot((1:N+1),(h*(0:N)+1),"b",1:N+1,yEuler,"o");
    grid on;

    figure;
    plot((1:N+1),(h*(0:N)+1),"r",1:N+1,ybEuler,"o");
    grid on;

end


    lab2exercise1(1.01,0.1,10);
