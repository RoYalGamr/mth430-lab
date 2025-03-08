clc;
clearvars;
%  We are solv1ng a 2 po1nt BVP problem w1th h3lp of Newtons method and RK4 method
% given problem
% u'' = -u + f   -1<t<1
% u(-1) = u(1) = 0
% where f is a some function of t , u , u'

% Initial S0
s2 = 0.2;

% Function f
f = @(t,y,dy) t + y + dy;
dfdy = @(t,y,dy) 1;
dfddy = @(t,y,dy) 1;

% Initial value according to the problem
a0=1;
a1=0;
b0=1;
b1=0;
c0=0;
c1=-1;
a=0;
b=0;
% y'' = m(y,y')
% z'' = l(y,y',z,z')
m = @(t,y,dy) ((-y)  + f(t,y,dy)) ;
l = @(t,y,dy,z,dz) ((-1 + dfdy(t,y,dy))*z + (dfddy(t,y,dy))*dz);

M = 12;
error = zeros(1:M);

for po = 3:M
    % Number of steps
    N = 2^po;
    h = 2/N;

    yexact = zeros(1,N+1);
    co1 = 2/(exp(1) - exp(-1));
    co2 = (-1/2) -co1*(exp(-1));
    for i = 1:N+1
        t = -1 + (i-1)*h;
        yexact(i) = -((t^2)/2) - t + (co1*exp(t)) + co2;
    end


    while (1)
        s1 = s2;

        y = zeros(N+1,4);

        % Initialised according to the problem
        y(1,:) = [a1*s1-c1*a, a0*s1-c0*a, a1, a0];

    %   General Code for every ODE
        % f = [y' , m , z' , l];
        F = @(t,y) [y(2), m(t,y(1),y(2)), y(4), l(t,y(1),y(2),y(3),y(4))];
        for i = 1:N
            t = -1 + (i-1)*h;
            k1 =  F(t, y(i,:));
            k2 = F(t+h/2, y(i,:) + h*k1/2);
            k3 = F(t+h/2, y(i,:) + h*k2/2);
            k4 = F(t+h, y(i,:) + h*k3);

            y(i+1,:) = y(i,:) + h*(k1 + (2*k2) + (2*k3) + k4)/6 ;
        end

    %  Dont forget to subtract a
        H = y(N+1,1) - b;
        dH = y(N+1,3);

        s2 = s1 - (H/dH);

        if (abs(s1-s2)<1e-8)
%             disp(s2);
%             t = linspace(-1, 1, N+1);
%             plot(t, y(:,1), '-');
%             xlabel('t');
%             ylabel('u(t)');
%             title('Solution of the BVP');
%             grid on;
            break;
        end
    end

    error(po) = max(abs(yexact-y(:,1)'))/max(abs(yexact));
end

fprintf('\n');
for m = 3:M
    fprintf('%6d \t %0.6e \t %0.2f\n', ...
    2^m, error(m), error(m-1)/error(m));
end

plot(-1+[0:N]*h,yexact,'LineWidth',2);
hold on;
plot(-1+ [0:N]*h,y(:,1)','--','LineWidth',2);
legend('yexact','y');
