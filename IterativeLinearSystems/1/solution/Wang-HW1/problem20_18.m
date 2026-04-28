A = [0.1, 0.5, -0.1; 0.4, 0.2, 0.6; 0.2, -0.3, 0.4];
sizeA = size(A, 1);

U = triu(A, 1);
D = diag(diag(A));
L = tril(A, -1);

BGS = -inv(L + D) * U;
rhoGS = max(abs(eig(BGS)));
n1_GS = norm(BGS, 1);
n2_GS = norm(BGS, 2);
nin_GS = norm(BGS, inf);
fprintf('1-norm of BGS = %.4f\n', n1_GS);
fprintf('2-norm of BGS = %.4f\n', n2_GS);
fprintf('infinity-norm of BGS = %.4f\n', nin_GS);
fprintf('Spectral radius of BGS = %.4f\n', rhoGS);

BJ = -inv(D) * (L + U);
n1_J = norm(BJ, 1);
n2_J = norm(BJ, 2);
nin_J = norm(BJ, inf);
rhoBJ = max(abs(eig(BJ)));
fprintf('1-norm of BJ = %.4f\n', n1_J);
fprintf('2-norm of BJ = %.4f\n', n2_J);
fprintf('infinity-norm of BJ = %.4f\n', nin_J);
fprintf('Spectral radius of BJ = %.4f\n', rhoBJ);

tol = 1.0e-14;
numiter = 2000000;
b = ones(sizeA, 1);
x0 = zeros(sizeA, 1);
[~, GSiter, GSresider] = gausseidel(A, b, x0, tol, numiter);
if GSiter > 0
  fprintf('Gauss-Seidel converges at step %.4f with residual %.4f\n', ...
          GSiter, GSresider);
else
  fprintf('Gauss-Seidel  does not converge\n');
end
[~, Jiter, Jresider] = Jacobi(A, b, x0, tol, numiter);
if Jiter > 0
  fprintf('Jacobi converges at step %.4f with residual %.4f\n', ...
          Jiter, Jresider);
else
  fprintf('Jacobi does not converge\n');
end
