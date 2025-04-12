clear all;

M = input(' ');
for m = 1:M
    N = 2^m; 

    h = 1/N;
    u = zeros(N-1,N-1);    % for storing the numerical solution
    ex = zeros(N-1,N-1);   % for storing the exact solution

    %%
    %% Complete the code to file u (numerical solution) and 
    %% ex (the exact solution)
    xv = linspace(0,1,N+1);
    yv = linspace(0,1,N+1);

    %% Exact Solution
    for j = 1:(N-1)
        for i = 1:(N-1)
            x = xv(i+1);
            y = yv(j+1);
            ex(j,i) = (((x^2) + (y^2))^2)/16;
        end
    end

    %% Numerical Solution

    nInterior = (N-1)^2;
    g = @(x,y) (((x.^2) + (y.^2)).^2)/16;

    A = sparse(nInterior, nInterior);
    b = zeros(nInterior, 1);

    for j = 1:(N-1)        
        for i = 1:(N-1)     
            k = (j-1)*(N-1) + i;

            x_ij = xv(i+1);
            y_ij = yv(j+1);

            % b(k) = x^2+y^2
            b(k) = (x_ij^2) + (y_ij^2);

            % Coefficient for u(i,j)
            A(k,k) = -4/(h^2);
    
            % Right neighbor (i+1,j)
            if i+1 <= (N-1)
                kRight = (j-1)*(N-1) + (i+1);
                A(k, kRight) = 1/(h^2);
            else
                % On boundary x = 1
                b(k) = b(k) - g(xv(end), y_ij) / (h^2);
            end
    
            % Left neighbor (i-1,j)
            if i-1 >= 1
                kLeft = (j-1)*(N-1) + (i-1);
                A(k, kLeft) = 1/(h^2);
            else
                % On boundary x = 0
                b(k) = b(k) - g(xv(1), y_ij) / (h^2);
            end
    
            % Top neighbor (i,j+1)
            if j+1 <= (N-1)
                kUp = (j)*(N-1) + i;
                A(k, kUp) = 1/(h^2);
            else
                % On boundary y = 1
                b(k) = b(k) - g(x_ij, yv(end)) / (h^2);
            end
    
            % Bottom neighbor (i,j-1)
            if j-1 >= 1
                kDown = (j-2)*(N-1) + i;
                A(k, kDown) = 1/(h^2);
            else
                % On boundary y = 0
                b(k) = b(k) - g(x_ij, yv(1)) / (h^2);
            end
   
        end
    end
    
    % Solve the linear system
    U_vec = A \ b;
    
    for j = 1:(N-1)
        for i = 1:(N-1)
            k = (j-1)*(N-1) + i;
            u(j,i) = U_vec(k);
        end
    end

    %%

    e(m) = max(max(abs(u-ex)));
end

for m = 1:M
    fprintf('%.2g ', e(m));    
end
fprintf('\n');

for m = 2:M
    fprintf('%.2g ', e(m-1)/e(m));    
end
fprintf('\n');
