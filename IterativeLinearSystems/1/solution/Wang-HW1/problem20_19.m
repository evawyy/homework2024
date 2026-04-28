n = 100;
r = 0.25;
a = (-r) * ones(n - 1, 1);  % 下对角
b =  (1 + 2 * r) * ones(n, 1);    % 主对角
c = (-r) * ones(n - 1, 1);  % 上对角

A = diag(b) + diag(a, -1) + diag(c, 1);
b = ones(n, 1);
tol = 0.5e-14;
numiter = 100;
x0 = zeros(n, 1);

[~, GSiter, GSresider] = gausseidel(A, b, x0, tol, numiter);
GSlogre = log(GSresider);
if GSiter > 0
  fprintf('Gauss-Seidel  converges at step %f with residual %f and log-residual %f\n', ...
          GSiter, GSresider, GSlogre);
else
  fprintf('Gauss-Seidel  does not converge\n');
end
[~, Jiter, Jresider] = Jacobi(A, b, x0, tol, numiter);
Jlogre = log(Jresider);
if Jiter > 0
  fprintf('Jacobi converges at step %f with residual %f and log-residual %f\n', ...
          Jiter, Jresider, Jlogre);

else
  fprintf('Jacobi does not converge\n');
end
W = [1.1, 1.2, 1.3, 1.5, 1.9];
for w = W
  [~, Soriter, Sorresider] = sor(A, b, x0, w, tol, numiter);
  Sorlogre = log(Sorresider);
  if Soriter > 0
    fprintf('SOR(%.1f) converges at step %f with residual %f and log-residual %f\n', ...
            w, Soriter, Sorresider, Sorlogre);
  else
    fprintf('SOR(%.1f) does not converge\n', w);
  end
end
x_dir = A \ b;
x_dir_relresid = norm(b - A * x_dir) / norm(b);
log_x_dir_residual = log(x_dir_relresid);
fprintf('B\\c residual is %f and log-residual %f\n', ...
        x_dir_relresid, log_x_dir_residual);
