clc;
clearvars;
%  We are solv1ng a 2 po1nt BVP problem w1th h3lp of finite difference  scheme
% given problem
% u'' = -u + 2(u')^2/u    -1<t<1
% u(-1) = u(1) = (e+e^-1)^-1

% Initial values according to the problem
a=-1;
b=1;
alpha = 1/(exp(1) + (1/exp(1)));
beta = 1/(exp(1) + (1/exp(1)));
% df/dy and df/d(y')
dfdy = @(t,y,y1) -1-(2*(y1)^2)/(y*y);
dfddy = @(t,y,y1) 4*y1/y;
% f(t,u,u') Consider it as a vector function
f = @(t,y,y1)  -y + (2*((y1).^2))./y;

M = 10;
error = zeros(1:M);

for po = 3:M
    % Number of steps
    N = 2^po;
    printf("Current N: %d",N);

    % Step size
    h = (b-a)/(N+1);

    % A matrix tridiagonal
    A = spdiags(ones(N,1)*[1 -2 1], -1:1, N, N);

    % Exact Solution
    yexact = zeros(1,N);
    for i = 2:N+1
        t = -1 + (i-1)*h;
        yexact(i-1) = (exp(t) + exp(-t))^-1;
    end

    % Initialise values;
    t = a + h*([1:N]');
    y = zeros(N,1);
    dy = zeros(N,1);
    ddy = zeros(N,1);
    g = zeros(N,1);
    g(1) = -alpha/(h*h);
    g(end) = -beta/(h*h);

    %  Initialized y with line interpolation
    for i = 1:N
        y(i) = alpha + i*(beta-alpha)/(N+1);
    end

    while (1)

        % y' Calculation
        dy(2:end-1) = (y(3:end) - y(1:end-2)) / (2*h);
        % Handle boundaries
        dy(1) = (y(2) - alpha) / (2*h);
        dy(end) = (beta - y(end-1)) / (2*h);

        % F vector
        fvec = f(t,y,dy);

        % We are solving Ay/(h*h) - fvec - g = 0
        H = (A*y)./(h*h) - fvec - g ;

        % Calculating F jacobian
        Fj = zeros(N,N);
        for i = 1:N
            Fj(i,i) = dfdy(a+i*h,y(i),dy(i));
        end
        for i = 1:N-1
            Fj(i,i+1) = dfddy(a+i*h,y(i),dy(i))/(2*h);
            Fj(i+1,i) = -dfddy(a+i*h,y(i+1),dy(i+1))/(2*h);
        end
        dH = A/(h*h) - Fj;

        % Newtons step
        y1 = y - (dH\H);

        % Epsilon check
        if (abs(norm(y1-y))<1e-8)
              y = y1;
              disp(y);
            break;
        end
        y = y1;
    end

    error(po) = max(abs(yexact-y'))/max(abs(yexact));
end

fprintf('\n');
for m = 3:M
    fprintf('%6d \t %0.6e \t %0.2f\n', 2^m, error(m), error(m-1)/error(m));
end

% Ploting
plot(a+[1:N]*h,yexact,'LineWidth',2);
hold on;
plot(a+[1:N]*h,y,'--','LineWidth',2);
legend('yexact','y');
