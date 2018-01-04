function [w, loss, yhat, xi] = svm_linf(X, y, C)
% Computes a soft margin classifier 
% Using primal problem
% assumes y a row vector in {-1,1}
% X is has each column a sample
%
% Optimizes:
%   min_w 1/n \sum_i hinge_loss(xi,yi;w) s.t. ||w||_inf <= C
%
%   Uses linear program to solve:
%       min 1/n sum_i -yi(w^T xi)   s.t. -C <= w_j <= C
%
%   Default: C = 1
%

if nargin ~= 3
    error('svm_linf requires exactly 3 arguments');
end

if C <= 0, error('C must be positive'); end

% if (nargin < 3 || C <= 0)
%     C = 1;
% end

X = X';             % now each row of X is an example
y = y';             % now y is a column vector in {-1, 1}
t = size(X,1);      % number of training examples
n = size(X,2);      % number of features 
% f = sum_{i=1}^t y_ix_i, assuming y_i a scalar
% f = -sum(X.*repmat(y',n,1)',1);
% f = f';

f = [zeros(n, 1); ones(t, 1)];

A = [repmat(-y, 1, n) .* X, -eye(t)];
b = -1 * ones(t, 1);

% linprog
lb = [-C * ones(n,1); zeros(t, 1)];
ub = [C * ones(n,1); inf * ones(t, 1)];
init_value = [zeros(n, 1); ones(t, 1)];
[w_xi,lploss,flag] = linprog(f, A, b, [], [], lb, ub, init_value, ...
                        optimset('Display', 'off', 'TolFun', 1e-4, 'TolX',1e-4));

loss = lploss;
% w is actually the first n entries in w
w = w_xi(1:n);
F = X * w;
yhat = sign(F);
xi = max(1 - y .* F, 0);

