clc;
clear;


function [E] = lab1exercise1(J)
#    y''' + 4y'' + 5y' + 2y = −4 sin t − 2 cost
#    y(0) = 1, y'(0) = 0, y''(0) = −1
    A = [0 1 0; 0 0 1; -2 -5 -4];

    E = zeros(J+1,1);

    for i = 0:J
        N = 2^i;
        Y = [1;0;-1];
        for n = 1:N
            h = 1/N;
            t = n/N;
            Y = Y + (h*(A*Y + [0;0;(-4*sin(t)-2*cos(t))]));
        end
        E(i+1) = abs(cos(1)- Y(1));
    end

end

E = lab1exercise1(10);
for j = 1:10
    disp(E(j)/E(j+1));
end

# The factor is approaching 2 .... But why ????


