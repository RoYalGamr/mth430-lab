T = 1.0;

M = 12;
error = zeros(1:M);

for m = 3:M
    N = 2^m;
    h = T/N;

    y = zeros(1,N+1);   % this array holds the exact solution
    yAB = zeros(1,N+1); % ... holds the numerical solution

    %%
    % Fill the arrays y and yAB appropriately

    % exact Solution
    for i = 1:N+1
        t = (i-1)*h;
        tmp = (5*(exp(-8*t/15)/6));
        y(i) = (1/3  - (tmp/5))/(1-tmp);
    end

    f = @(y) -4*(y-(1/5))*(y-(1/3));
    ff = @(y) -4*(2*y-(8/15))*(f(y));
    fff = @(y) 16*f(y)*((2*(y-1/5)*(y-1/3)) + ((2*y-8/15)*(y-1/3)) + ((2*y-8/15)*(y-1/5)));

    yAB(1) = 1;
    for i = 1:3
        yAB(i+1) = yAB(i) + (h*f(yAB(i))) + (h*h*ff(yAB(i))/2) + (h*h*h*f(yAB(i))/6);
    end

    for i = 4:N
        b0 = 55*(f(yAB(i)))/24;
        b1 = -59*(f(yAB(i-1)))/24;
        b2 = 37*(f(yAB(i-2)))/24;
        b3 = -9*(f(yAB(i-3)))/24;
        yAB(i+1) = yAB(i) + h*(b0+b1+b2+b3);
    end

    % Do not change the template below this line
    %%


    error(m) = max(abs(y-yAB))/max(abs(y));
end

fprintf('\n');
for m = 3:M
    fprintf('%6d \t %0.6e \t %0.2f\n', ...
    2^m, error(m), error(m-1)/error(m));
end

plot([0:N]*h,y,'LineWidth',2);
hold on;
plot([0:N]*h,yAB,'--','LineWidth',2);
legend('y','yAB');
