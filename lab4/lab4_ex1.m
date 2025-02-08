function lab4_ex1(y0,h,N)
#y' = −100y + 100t + 101,
#y(0) = y0
    y = zeros(N+1,1);
    yM1 = zeros(N+1,1);
    yM2 = zeros(N+1,1);

    yM1(1) = yM2(1) = y0;

    % Exact Solution
    for i = 0:N
        y(i+1) = i*h + 1 + (y0-1)*(exp(-100*i*h));
    end

    % Here comes our lamda functions
    f = @(t,y) -100*y+100*t+101;
    ft = @(t,y) 100;
    fy = @(t,y) -100;

    % M3th0d 1
    for i = 1:N
        t = i*h;
        yM1(i+1) = yM1(i) + (h*f(t,yM1(i))) + ((h^2)*(ft(t,yM1(i)) + fy(t,yM1(i))*f(t,yM1(i))))/2;
    end

    ftt = @(t,y) 0;
    fyy = @(t,y) 0;
    fty = @(t,y) 0;

    % M3th0d 2
    for i = 1:N
        t = i*h;
        yM2(i+1) = yM2(i) + h*(f(t,yM2(i))) + ((h^2)*(ft(t,yM2(i)) + fy(t,yM2(i))*f(t,yM2(i))))/2;
        yM2(i+1) = yM2(i+1) + ((h^3)*(ftt(t,yM2(i)) + 2*fty(t,yM2(i)) + fyy(t,yM2(i))*f(t,yM2(i))*f(t,yM2(i)) +ft(t,yM2(i))*fy(t,yM2(i)) + fy(t,yM2(i))*fy(t,yM2(i))*f(t,yM2(i)) ))/6;
    end

    % Plot
    figure;
    t = h*[0:N];
    plot(t,y,'o',t,yM1,'-',t,yM2,"g");
    legend("exact","Method 1","Method 2");

    % Error Calculations
    E1 = abs(yM1 - y);
    E2 = abs(yM2 - y);
    printf("%e %e\n",E1,E2);
end
