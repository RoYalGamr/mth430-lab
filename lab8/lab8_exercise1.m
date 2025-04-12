clc; clear; close all;

% Domain and grid parameters
N = 128;               % N+1 grid points in each coordinate (0 to 1)
h = 2.0 / N;          % grid spacing
xv = linspace(-1,1,N+1);
yv = linspace(-1,1,N+1);
nInterior = (N-1)^2;  % number of unknowns (only interior points)

% Dirichlet boundary function: g(x,y)
g = @(x,y) log((x-3).^2 + (y-2).^2);

% Assembling the linear system: A*U_vec = b.
% U_vec holds u at interior points in row-major order.
A = sparse(nInterior, nInterior);
b = zeros(nInterior, 1);

for j = 1:(N-1)         % interior y-index (j=1 corresponds to y=xv(2))
    for i = 1:(N-1)     % interior x-index (i=1 corresponds to xv(2))
        k = (j-1)*(N-1) + i;
        % Coordinates at interior point (i,j)
        x_ij = xv(i+1);
        y_ij = yv(j+1);

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

        % Right-hand side f(x,y)=0 in the interior for Laplace's eqn.
        % b(k) remains unchanged.
    end
end

% Solve the linear system
U_vec = A \ b;

% Map solution to full grid (include boundary)
U_mat = zeros(N+1, N+1);
for j = 1:(N-1)
    for i = 1:(N-1)
        k = (j-1)*(N-1) + i;
        U_mat(j+1,i+1) = U_vec(k);
    end
end

% Set boundary values using g(x,y)
U_mat(1,:)   = g(xv, yv(1));       % bottom (y=0)
U_mat(end,:) = g(xv, yv(end));       % top (y=1)
U_mat(:,1)   = g(xv(1), yv);         % left (x=0)
U_mat(:,end) = g(xv(end), yv);         % right (x=1)

[X, Y] = meshgrid(xv, yv);
figure; surf(X, Y, U_mat); shading interp; colorbar;
xlabel('x'); ylabel('y'); zlabel('u(x,y)');
title('Solution with u = log((x-3)^2+(y-2)^2) on the boundary');

