clc;
clear;
clearvars;
function [y,yEX,yEM] = lab3_ex1(y0,h,N)
    y = zeros(N+1,1);
    yEX = zeros(N+1,1);
    yEM = zeros(N+1,1);


    yEM(1) = yEX(1) = y0;

    % Exact Solution
    for i = 0:N
        y(i+1) = i*h + 1 + (y0-1)*(exp(-100*i*h));
    end

    yEM(2) = y0 + h*(-100*y0 + 101);
    for i = 2:N
        fnn = ((-100*yEM(i-1)) + (100*(i-2)*h) +101);
        fn = ((-100*yEM(i)) + (100*(i-1)*h) +101);
        yEM(i+1) = yEM(i) + (2*fn*h/3) - (fnn*h/12) + (5*h*(100*((i-1)*h)+101)/12);
        yEM(i+1) = yEM(i+1)/(1+(125*h/3));
    end


    yEX(2) = y(2);
    for i = 2:N
        fnn = ((-100*yEX(i-1)) + (100*(i-2)*h) +101);
        fn = ((-100*yEX(i)) + (100*(i-1)*h) +101);
        yEX(i+1) = yEX(i) + (2*fn*h/3) - (fnn*h/12) + (5*h*(100*((i-1)*h)+101)/12);
        yEX(i+1) = yEX(i+1)/(1+(125*h/3));
    end
    t = (0:N) * h;
    plot(t,y,'g',t,yEX,'r',t,yEM,'b');
    legend('Exact Solution','EM','EX');

end


[y,yEX,yEM] = lab3_ex1(-1,0.0001,500);
