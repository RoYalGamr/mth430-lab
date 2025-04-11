clc; clear; close all;
N = 128;
h = 1.0/N;
xv = linspace(0,1,N+1);
yv = linspace(0,1,N+1);
g = @(x,y) log((x-3).^2 + (y-2).^2);
inDomain = @(x,y) ~((x<=0.5) & (y<=0.5));
indexInterior = zeros(N-1, N-1);
idx = 0;
for j = 1:(N-1)
    for i = 1:(N-1)
        x_val = xv(i+1);
        y_val = yv(j+1);
        if inDomain(x_val,y_val)
            idx = idx + 1;
            indexInterior(j,i) = idx;
        end
    end
end
nInterior = idx;
A = sparse(nInterior, nInterior);
b = zeros(nInterior,1);
for j = 1:(N-1)
    for i = 1:(N-1)
        k = indexInterior(j,i);
        if k==0, continue; end
        x_val = xv(i+1);
        y_val = yv(j+1);
        A(k,k) = -4/(h^2);
        if i+1 <= (N-1)
            kR = indexInterior(j,i+1);
            if kR > 0
                A(k,kR) = 1/(h^2);
            else
                b(k) = b(k) - g(xv(i+2),y_val)/(h^2);
            end
        else
            b(k) = b(k) - g(xv(end),y_val)/(h^2);
        end
        if i-1 >= 1
            kL = indexInterior(j,i-1);
            if kL > 0
                A(k,kL) = 1/(h^2);
            else
                b(k) = b(k) - g(xv(i),y_val)/(h^2);
            end
        else
            b(k) = b(k) - g(xv(1),y_val)/(h^2);
        end
        if j+1 <= (N-1)
            kU = indexInterior(j+1,i);
            if kU > 0
                A(k,kU) = 1/(h^2);
            else
                b(k) = b(k) - g(x_val,yv(j+2))/(h^2);
            end
        else
            b(k) = b(k) - g(x_val,yv(end))/(h^2);
        end
        if j-1 >= 1
            kD = indexInterior(j-1,i);
            if kD > 0
                A(k,kD) = 1/(h^2);
            else
                b(k) = b(k) - g(x_val,yv(j))/(h^2);
            end
        else
            b(k) = b(k) - g(x_val,yv(1))/(h^2);
        end
    end
end
U_vec = A\b;
U_mat = nan(N+1, N+1);
for j = 1:(N-1)
    for i = 1:(N-1)
        k = indexInterior(j,i);
        if k > 0
            U_mat(j+1,i+1) = U_vec(k);
        end
    end
end
for i = 1:N+1
    if inDomain(xv(i),yv(1))
        U_mat(1,i) = g(xv(i),yv(1));
    end
    if inDomain(xv(i),yv(end))
        U_mat(end,i) = g(xv(i),yv(end));
    end
end
for j = 1:N+1
    if inDomain(xv(1),yv(j))
        U_mat(j,1) = g(xv(1),yv(j));
    end
    if inDomain(xv(end),yv(j))
        U_mat(j,end) = g(xv(end),yv(j));
    end
end
for j = 1:N+1
    for i = 1:N+1
        if inDomain(xv(i),yv(j)) && isnan(U_mat(j,i))
            U_mat(j,i) = g(xv(i),yv(j));
        end
    end
end
[X,Y] = meshgrid(xv,yv);
figure; surf(X,Y,U_mat); shading interp; colorbar;
xlabel('x'); ylabel('y'); zlabel('u(x,y)');
title('Solution with u = log((x-3)^2+(y-2)^2)');

