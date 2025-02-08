function [Px, Py, Pz] = lab4_ex2(lat0, lon0, v0, gam0, phi0, h, T)
    % Constants
    r0 = 6378137;
    u = 3986004.418 * (10^8);
    j2 = 1091.3 * (10^-6);

    % Lambda functions
    fu = @(r, p, c) (-u * p / (r^3)) * ((r0 / r)^2 + (3/2) * j2 * ((r0 / r)^4) * (c - 5 * ((p / r)^2)));
    g = @(P) [fu(norm(P), P(1), 1), fu(norm(P), P(2), 1), fu(norm(P), P(3), 3)];
    F = @(P, V) [V, g(P)];

    % Initial Values
    P0 = [r0 * cos(lat0) * cos(lon0), r0 * cos(lat0) * sin(lon0), r0 * sin(lat0)];
    A = [-sin(lat0) * cos(lon0), sin(lon0), cos(lat0) * cos(lon0);
         -sin(lat0) * sin(lon0), -cos(lon0), cos(lat0) * sin(lon0);
          cos(lat0), 0, sin(lat0)];
    V0 = A * [v0 * cos(gam0) * cos(phi0); v0 * cos(gam0) * sin(phi0); v0 * sin(gam0)];

    % RK4 Method
    N = T / h;
    P = zeros(N+1, 3);
    V = zeros(N+1, 3);
    P(1, :) = P0;
    V(1, :) = V0';

    for i = 1:N
        k1 = F(P(i, :), V(i, :));
        k2 = F(P(i, :) + h/2 .* k1(1:3), V(i, :) + h/2 .* k1(4:6));
        k3 = F(P(i, :) + h/2 .* k2(1:3), V(i, :) + h/2 .* k2(4:6));
        k4 = F(P(i, :) + h .* k3(1:3), V(i, :) + h .* k3(4:6));

        P(i+1, :) = P(i, :) + h/6 .* (k1(1:3) + 2*k2(1:3) + 2*k3(1:3) + k4(1:3));
        V(i+1, :) = V(i, :) + h/6 .* (k1(4:6) + 2*k2(4:6) + 2*k3(4:6) + k4(4:6));
    end

    % Plotting Sphere and Trajectory
    figure;
    theta = linspace(0, 2*pi, 100);
    phi = linspace(0, pi, 100);
    [thetaGrid, phiGrid] = meshgrid(theta, phi);

    xSphere = r0 * sin(phiGrid) .* cos(thetaGrid);
    ySphere = r0 * sin(phiGrid) .* sin(thetaGrid);
    zSphere = r0 * cos(phiGrid);

    Px = P(:, 1);
    Py = P(:, 2);
    Pz = P(:, 3);

    surf(xSphere, ySphere, zSphere, 'EdgeColor', 'none', 'FaceAlpha', 0.7); % Sphere
    hold on;
    plot3(Px, Py, Pz, 'r', 'LineWidth', 2); % Trajectory
    xlabel('X Coordinate');
    ylabel('Y Coordinate');
    zlabel('Z Coordinate');
    title('Trajectory on Sphere');
    axis equal;
    grid on;
    rotate3d on;
end

