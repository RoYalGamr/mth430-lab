clear all;

k = -14;
M = 12;

for m = 1:M
    n = 2^m;
    h = 2/(n+1);

    u = zeros(n,1);     % for storing the numerical solution
    
    %%
    %% complete the code to fill the vector u
    f = @(t,y,dy) ((t*dy)-(k*k*y))/(1-(t*t));
    f2 = @(t,y,dy) (-(k*k)/(1-(t*t)));
    f3 = @(t,y,dy) (t/(1-(t*t)));
    f_vec = zeros(n,1);
    du = zeros(n,1);
    alpha = (-1)^k;
    beta = 1;
    a = -1;
    b = 1;

    % Initialise u with line interpolation
    for i = 1:n
        u(i) = alpha + i*(beta-alpha)/(n+1);
    end

    % matrix A creation
    A = zeros(n,n);
    for i = 1:n
        A(i,i) = -2;
    end
    for i = 1:n-1
        A(i+1,i) = 1;
        A(i,i+1) = 1;
    end
    
    % g vector creation
    g = zeros(n,1);
    g(1) = -alpha/(h*h);
    g(end) = -beta/(h*h);


    while (1)
        % dy calculation
        for i = 2:n-1
            du(i) = (u(i+1) - u(i-1))/(2*h);
        end
        du(1) = (u(2) - alpha)/(2*h);
        du(n) = (beta - u(n-1))/(2*h);
        
        % f vec calculation
        for i = 1:n
            f_vec(i) = f(a+i*h,u(i),du(i));
        end

        H = ((A*u)./(h*h)) - f_vec - g;
        
        % f vec jacobian calculation
        Fj = zeros(n,n);
        for i = 1:n
            Fj(i,i) = f2(a+i*h,u(i),du(i));
        end
        for i = 1:n-1
            Fj(i+1,i) = -f3(a+(i+1)*h,u(i+1),du(i+1))/(2*h);
            Fj(i,i+1) = f3(a+i*h,u(i),du(i))/(2*h);
        end

        dH = (A./(h*h))  - Fj;
        
        % Newton step
        u1 = u - (dH\H);

        % Epsilon check
        if (abs(norm(u1 - u))< 1e-10)
            u = u1;
            break;
        end
        u = u1;
    end
    
    %%
    
    ue = cos(k*acos([-1+h:h:1-h].'));   % exact solution u
    E(m) = max(abs(u-ue));              % error
end

fprintf('\n');
for m = 3:M
    fprintf('%6d \t %0.6e \t %0.2f\n', ...
    2^m, E(m), E(m-1)/E(m));
end

plot([0:n+1]*h-1,[(-1)^k;ue;1],'LineWidth',2);
hold on;
plot([0:n+1]*h-1,[(-1)^k;u;1],'--','LineWidth',2);
legend('ue','u');
