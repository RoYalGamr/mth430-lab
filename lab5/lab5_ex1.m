clc;
clearvars;
%  We are solv1ng a 2 po1nt BVP problem w1th h3lp of Newtons method and RK4 method
% given problem
% u'' = -u + f   -1<t<1
% u(-1) = u(1) = 0
% where f is a some function of t , u , u'

% Number of steps
N = 10000;
h = 2/N;

% Initial S0
s2 = 0.2;

% Function f
f = @(t,y,dy) t + y + dy;
dfdy = @(t,y,dy) 1;
dfddy = @(t,y,dy) 1;

a =  0;

while (1)
    s1 = s2;

    y = zeros(N+1,4);

    % Initialised according to the problem
    y(1,:) = [a, s1, 0, 1];
    % y'' = m(y,y')
    % z'' = l(y,y',z,z')
    m = @(t,y,dy) ((-y)  + f(t,y,dy)) ;
    l = @(t,y,dy,z,dz) ((-1 + dfdy(t,y,dy))*z + (dfddy(t,y,dy))*dz);

%   General Code for every ODE
    % f = [y' , m , z' , l];
    F = @(t,y) [y(2), m(t,y(1),y(2)), y(4), l(t,y(1),y(2),y(3),y(4))];
    for i = 1:N
        t = (i-1)*h;
        k1 =  F(t, y(i,:));
        k2 = F(t+h/2, y(i,:) + h*k1/2);
        k3 = F(t+h/2, y(i,:) + h*k2/2);
        k4 = F(t+h, y(i,:) + h*k3);

        y(i+1,:) = y(i,:) + h*(k1 + (2*k2) + (2*k3) + k4)/6 ;
    end

%  Dont forget to subtract a
    H = y(N+1,1) - a;
    dH = y(N+1,3);


    s2 = s1 - (H/dH);

    if (abs(s1-s2)<1e-8)
        disp(s2);
        t = linspace(-1, 1, N+1);
        plot(t, y(:,1), '-');
        xlabel('t');
        ylabel('u(t)');
        title('Solution of the BVP');
        grid on;
        break;
    end
end


