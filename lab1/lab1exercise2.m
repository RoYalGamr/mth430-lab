clc;
clearvars;
clear;

    function [] = lab1exercise2(N, T , a1, a2, b1, b2, y10, y20)
        E = zeros(N+1,2);
        E(1,1) = y10;
        E(1,2) = y20;
        A = [a1 b2;-b1 -a2];

        for n = 1:N
            t = n*T/N;
            h = T/N;
            E(n+1,1) = E(n,1) + h*(E(n,1)*(a1 - b1*E(n,2)));
            E(n+1,2) = E(n,2) + h*(E(n,2)*(-a2 +b2*E(n,1)));
        end
        disp(E);

        figure;
        plot(1:N+1,E(:,1),"r",1:N+1,E(:,2),"b");
        xlabel("time");
        ylabel("population");

        figure;
        plot3(1:N+1,E(:,1),E(:,2));
        xlabel("time");
        ylabel("prey");
        zlabel("predator");
        rotate3d;
        grid on;
    end

lab1exercise2(1000,25,1,0.5,0.1,0.02,100,10);
