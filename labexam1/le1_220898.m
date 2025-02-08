T = 1;

M = 10;
for m = 1:M
    N = 2^m;
    h = T/N;

    y = zeros(1,N+1);  % this array holds the exact solution
    yE = zeros(1,N+1); % holds the numerical solution obtained using Euler's method
    yB = zeros(1,N+1); % holds the numerical solution obtained using Backward Euler
    yT = zeros(1,N+1); % holds the numerical solution obtained using Trapezoidal method

    %%
    % Fill the arrays y, yE, yB and yT appropriately


    % Exact solution
    for t = 0:N
        x = t*h;
        y(t+1) = exp(2*sin(20*x)) + exp(-50*x);
    end

    % Euler method
    yE(1) = 2;

    for t = 1:N
        x = (t-1)*h;
        yE(t+1) = yE(t) + h*(-50*(yE(t) - exp(2*(sin(20*x)))) + 40*cos(20*x)*exp(2*(sin(20*x))));
    end

    % Backward Euler method
    yB(1) = 2;

    for t = 1:N
        x = t*h;
        yB(t+1) = yB(t) + h*(50*(exp(2*(sin(20*x)))) + 40*cos(20*(x))*exp(2*(sin(20*x))));
        yB(t+1) = yB(t+1)/(1+50*h);

    end

    % Trapezoid Method
    yT(1) = 2;

    for t = 1:N
        x1 = (t-1)*h;
        x2 = t*h;
        yT(t+1) = yT(t) + h*(-50*(yT(t) - exp(2*(sin(20*x1)))) + 40*cos(20*x1)*exp(2*(sin(20*x1))))/2;
        yT(t+1) = yT(t+1) + h*(50*(exp(2*(sin(20*x2)))) + 40*cos(20*(x2))*exp(2*(sin(20*x2))))/2;
        yT(t+1) = yT(t+1)/(1+(50*h)/2);

    end


    % Do not change the template below this line
    %%

    eE(m) = max(abs(yE-y))/max(abs(y));
    eB(m) = max(abs(yB-y))/max(abs(y));
    eT(m) = max(abs(yT-y))/max(abs(y));
end

% Convergence report
for m = 2:M
    fprintf('%4d \t %.2e \t %.1e \t %.2e \t %.1e \t %.2e \t %.1e \n', ...
        2^m, eE(m), eE(m-1)/eE(m), eB(m), eB(m-1)/eB(m), ...
        eT(m), eT(m-1)/eT(m));
end

#plot([0:N]*h,y,'LineWidth',2);
hold on;
#plot([0:N]*h,yE,'--','LineWidth',2);
#plot([0:N]*h,yB,'-.','LineWidth',2);
plot([0:N]*h,yT,':','LineWidth',2);
legend('y','yE','yBE','yT');
