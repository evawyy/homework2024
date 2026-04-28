size=100;
b=ones(size,1);
A=diag(1:size)+diag(ones(size-1,1),1)+diag(ones(size-1,1),-1);
x0=zeros(size,1);
tol=1e-14;
maxiter=100;

kappa=cond(A);
q=(sqrt(kappa)-1)/(sqrt(kappa)+1);
 function [x, residual, iter,rs,true_res,bound] = cg(A,b,x0,tol,maxiter,q)
%cg  Conjugate gradient method
%
%   [x,residual,iter] = cg(A,b,x0,tol,numitr) computes the solution x 
%   of a sparse symmetric positive definite system Ax = b using
%   the Conjugate Gradient method. x0 is the initial approximation, tol is
%   the error tolerance, and maxiter is the maximum number of iterations
%   to execute. If the conjugate gradient method converges, iter contains
%   the number of iterations required to converge; otherwise, iter = -1.
%   For consistent failure of convergence, try precg.
rs=zeros(maxiter+1,1);
true_res=zeros(maxiter+1,1);
bound=zeros(maxiter+1,1);
r = b - A*x0;
p = r;
x = x0;
tol = tol^2;
normr = norm(r);
rs(1)=normr;
true_res(1)=norm(b-A*x);
bound(1)=2;
for i = 1:maxiter
  normrsqr = r'*r;
   w = A*p;
   alpha = normrsqr/(p'*w);
   x = x + alpha*p;
   r = r - alpha*w;
   normrnewsqr = r'*r;
   if normrnewsqr < tol
      iter = i;
      residual = sqrt(normrnewsqr);
      return;
   end
   beta = normrnewsqr/normrsqr;
   p = r + beta*p;
    normr=norm(r);
   rs(i+1)=normr;
   true_res(i+1)=norm(b-A*x);
    bound(i+1)=2*(q^i);
end

iter = -1;
residual = sqrt(normrnewsqr);
end

function [x, r, iter,rs] = steepestDescent(A,b,x0,tol,maxiter)
%steepestDescent Steepest descent method for solving the sparse
%positive definite system Ax = b.
%
%   [x,iter] = steepestDescent(A,x0,b,tol,numitr) computes the solution x 
%   of a sparse symmetric positive definite system Ax = b using
%   the method of steepst descent. x0 is the initial approximation, tol is
%   the error tolerance, and maxiter is the maximum number of iterations
%   to execute. If the method converges, iter contains
%   the number of iterations required to converge; otherwise, iter = -1.
%   This function serves only as motivation for the conjugate gradient
%   method.

r = b - A*x0;
x = x0;
iter = 1;
rs=zeros(maxiter,1);
rs(1)=r'*r;
while (norm(r) >= tol) && (iter <= maxiter)

   alpha = (r'*r)/(r'*(A*r));
   x = x + alpha*r;
   r = b - A*x;
   iter = iter + 1;

    normr=norm(r);
    rs(iter)=normr;
end
iter = iter-1;

if iter >= maxiter
   iter = -1;
end
end 
[~, sd_r, sd_iter,sd_res] = steepestDescent(A,b,x0,tol,maxiter);
[~, cg_residual, cg_iter,cg_res,cg_true,bound] = cg(A,b,x0,tol,maxiter,q);
nvals = 0:maxiter;
figure;

semilogy(nvals, cg_res, 'r', 'LineWidth', 1.5); hold on;
semilogy(nvals, cg_true, 'b--', 'LineWidth', 1.5);
semilogy(nvals, sd_res, 'g', 'LineWidth', 1.5);
semilogy(nvals, bound, 'k:', 'LineWidth', 1.5);

legend('CG computed','CG true','Steepest descent','CG bound');
xlabel('Iteration n');
ylabel('Residual norm');
title('Convergence behavior');
grid on;

exportgraphics(gcf, 'convergence.pdf');   
