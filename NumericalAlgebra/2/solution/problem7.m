% problem7_cgne.m
% Reproduce CGNE residual behavior figure

clear; clc; close all;

rng(0);

% Parameters
m = 40;
kappa = 1e4;
maxit = 600;

% Generate random orthogonal matrices U and V
[U, ~] = qr(randn(m));
[V, ~] = qr(randn(m));

% Singular values
sigma = kappa.^(-(0:m-1)/(m-1));

% Construct matrix A = U * Sigma * V'
A = U * diag(sigma) * V';

% Exact solution and right-hand side
xstar = randn(m, 1);
b = A * xstar;

% Initial guess
x = zeros(m, 1);

% CGNE initialization
r = b - A * x;
z = A' * r;
p = z;

% Storage
true_res = zeros(maxit, 1);
upd_res = zeros(maxit, 1);
xnorm_ratio = zeros(maxit, 1);

scale = norm(A) * norm(xstar);

% CGNE iteration
for k = 1:maxit
    Ap = A * p;

    alpha = (z' * z) / (Ap' * Ap);

    x = x + alpha * p;

    % Updated residual
    r = r - alpha * Ap;

    znew = A' * r;

    beta = (znew' * znew) / (z' * z);

    p = znew + beta * p;
    z = znew;

    % True residual recomputed directly
    true_res(k) = norm(b - A * x) / scale;

    % Updated residual
    upd_res(k) = norm(r) / scale;

    % Norm ratio
    xnorm_ratio(k) = norm(x) / norm(xstar);
end

% Find index where ||x_n|| / ||x*|| is maximal
[~, imax] = max(xnorm_ratio);

% Plot
figure;
semilogy(1:maxit, true_res, '-', 'LineWidth', 1.5);
hold on;
semilogy(1:maxit, upd_res, '--', 'LineWidth', 1.5);
semilogy(imax, true_res(imax), 'ro', 'MarkerSize', 8, 'LineWidth', 2);
grid on;

legend('true residual', ...
       'updated residual', ...
       'max ||x_n||/||x^*||', ...
       'Location', 'best');

xlabel('iteration n');
ylabel('relative residual');
title('CGNE residual behavior');

% Save figures
exportgraphics(gcf, 'cgne_residual.pdf', 'ContentType', 'vector');
saveas(gcf, 'cgne_residual.png');

fprintf('Figure saved as cgne_residual.pdf and cgne_residual.png\n');
fprintf('Maximum ||x_n||/||x^*|| occurs at iteration %d\n', imax);
fprintf('Maximum ratio = %.6e\n', xnorm_ratio(imax));
