clc;
clearvars;

% Problem setup
a = -1;
b = 1;
alpha = 1/(exp(1) + exp(-1));
beta = alpha;

% Functions
f = @(t,y,y1) -y + (2*(y1.^2))./y;
dfdy = @(t,y,y1) -1 - (2*(y1.^2))./(y.^2);
dfddy = @(t,y,y1) (4*y1)./y;

M = 10;
error = zeros(1,M);

% Use po = 3:M to compute errors for multiple N, or 4:4 for just N=16
for po = 3:8
    N = 2^po;
    fprintf("Current N: %d\n", N);
    h = (b-a)/N;
    t = a + h*(0:N)'; % Column vector

    % Precompute matrices
    V = t.^(0:N); % Vandermonde: V(i,j) = t(i)^(j-1)
    D = zeros(N+1, N+1);
    for j = 2:N+1
        D(:,j) = (j-1) * t.^(j-2);
    end

    % Matrix A
    A = zeros(N+1, N+1);
    A(1,:) = a.^(0:N); % u(-1)
    A(N+1,:) = b.^(0:N); % u(1)
    for i = 2:N
        for j = 3:N+1
            A(i,j) = (j-1)*(j-2)*(t(i)^(j-3));
        end
    end

    % Exact solution
    yexact = (exp(t) + exp(-t)).^-1;

    % Initial guess
%    y = (-alpha)* ones(N+1,1); % Matches boundary conditions
    tmp = [5.0000e-01,-1.5266e-12,-2.5000e-01, 1.2725e-12, 1.0417e-01,-8.0158e-13,-4.2361e-02, 1.3161e-12, 1.7175e-02,-1.7244e-11,-6.9611e-03, 2.0115e-10, 2.8212e-03,-1.5835e-09,-1.1434e-03, 8.6784e-09, 4.6339e-04,-3.3935e-08,-1.8776e-04, 9.5831e-08, 7.5963e-05,-1.9559e-07,-3.0491e-05, 2.8531e-07, 1.1879e-05,-2.8949e-07,-4.2565e-06, 1.9373e-07, 1.2689e-06,-7.6735e-08,-2.6684e-07, 1.3604e-08, 2.8428e-08];
    y = zeros(N+1,1);
    y(1:min(N+1,33)) = tmp'(1:min(N+1,33));
    % Newton's method
    iter = 0;
    max_iter = 100;
    while iter < max_iter
        % Compute u and u' efficiently
        sum1 = V * y;
        sum2 = D * y;

        % Residual vector
        fvec = zeros(N+1,1);
        fvec(1) = alpha;
        fvec(N+1) = beta;
        fvec(2:N) = f(t(2:N), sum1(2:N), sum2(2:N));
        H = A*y - fvec;

        % Jacobian Fj
        Fj = zeros(N+1, N+1);
        dfdy_vec = dfdy(t(2:N), sum1(2:N), sum2(2:N));
        dfddy_vec = dfddy(t(2:N), sum1(2:N), sum2(2:N));
        Fj(2:N,:) = dfdy_vec .* V(2:N,:) + dfddy_vec .* D(2:N,:);
        dH = A - Fj;

        % Check condition number
        cond_number = cond(dH);
        fprintf('Iteration %d, Condition number: %e\n', iter, cond_number);

        % ----- Adaptive lambda (line search) -----
        % Choose candidate lambda values (for instance, 50 values linearly spaced in [0.01, 1])
        candidateLambdas = linspace(0.001, 1, 500);
        bestLambda = candidateLambdas(1);
        bestNorm = inf;

        delta = dH \ H;
        % Test each candidate
        for L = candidateLambdas
            y_candidate = y - L * delta;

            % Compute residual for candidate:
            sum1_candidate = V * y_candidate;
            sum2_candidate = D * y_candidate;
            fvec_candidate = zeros(N+1,1);
            fvec_candidate(1) = alpha;
            fvec_candidate(N+1) = beta;
            fvec_candidate(2:N) = f(t(2:N), sum1_candidate(2:N), sum2_candidate(2:N));
            H_candidate = A*y_candidate - fvec_candidate;

            candidateNorm = norm(H_candidate);

            if candidateNorm < bestNorm
                bestNorm = candidateNorm;
                bestLambda = L;
            end
        end
        % Use the best lambda found
        lamda = bestLambda;
        fprintf('Iteration %d, selected lambda: %f, Residual norm: %e\n', iter, lamda, bestNorm);


        y1 = y - (lamda)*(delta);
        if norm(y1 - y) < 1e-6
            y = y1;

            break;
        end
        y = y1;
        iter = iter + 1;
    end
    if iter == max_iter
        warning('Maximum iterations reached');
    end
    disp(norm(H));
    y

    % Solution
    sol = V * y;
    error(po) = max(abs(yexact - sol))/max(abs(yexact));
end

% Display errors (adjust if po loop changes)
fprintf('\nN\tError\t\tRatio\n');
for m = 3:M
    if m <= po
        fprintf('%6d\t%0.6e\t%0.2f\n', 2^m, error(m), ifelse(m > 3, error(m-1)/error(m), NaN));
    end
end

% Plot
plot(t, yexact, 'LineWidth', 2);
hold on;
plot(t, sol, '--', 'LineWidth', 2);
legend('Exact', 'Computed');
