function [sol, loss] = svm_leastsq(X, y, C)

% X: #basis * #labeled
% y: 2 * #labeled

num_basis = size(X, 1);
U0 = zeros(num_basis, 1);
H = X * X';
f = - X * y(1, :)';

lb = -C * ones(num_basis, 1);
ub = C * ones(num_basis, 1);

[sol, fval, flag, iter] = quadprog(H, f, [], [], [], [], lb, ub, U0, ...
                            optimset('Display', 'off', 'TolFun', 1e-6, 'TolX',1e-6));

loss = fval * 2 + norm(y(1, :))^2;

sol = [sol -sol]';

end
