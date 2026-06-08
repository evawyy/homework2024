%a
function I = glquad(f,n)
beta = (1:n-1)./sqrt(4*(1:n-1).^2-1);
J = diag(beta,1) + diag(beta,-1);
[V,D] = eig(J);
x = diag(D); [x,idx] = sort(x); V = V(:,idx);
w = 2*(V(1,:).^2)';
I = w.'*f(x);
end
% -------- (b) --------
exact = exp(1)-exp(-1);
err = zeros(40,1);
for n = 1:40
    err(n) = abs(exact - glquad(@(x) exp(x), n));
end

figure;
semilogy(1:40, err, 'o-'); grid on;
xlabel('n'); ylabel('|I(e^x)-I_n(e^x)|');
title('Gauss-Legendre error for exp(x)');

exportgraphics(gcf, 'exp_error.pdf', 'ContentType','vector');

% -------- (c) --------
exact = 2*(exp(1)-1);
err = zeros(40,1);
for n = 1:40
    err(n) = abs(exact - glquad(@(x) exp(abs(x)), n));
end

figure;
semilogy(1:40, err, 'o-'); grid on;
xlabel('n'); ylabel('|I(e^{|x|})-I_n(e^{|x|})|');
title('Gauss-Legendre error for exp(abs(x))');

exportgraphics(gcf, 'exp_abs_error.pdf', 'ContentType','vector');
