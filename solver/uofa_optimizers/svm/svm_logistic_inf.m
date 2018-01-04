function [sol, loss] = svm_logistic_inf(X, y, C, verbose)

% X: #basis * #labeled
% y: 1 * #labeled

num_basis = size(X, 1);
num_task = size(y, 1);

y = y';

lb = -C * ones(num_basis, 1);
ub = C * ones(num_basis, 1);

w0 = zeros(num_basis, 1);

% [sol, loss, exitflag, output] = fmincon(@local_logistic_loss, w0, ...
%                                     [], [], [], [], lb, ub, [], ...
%                                     optimset('Display', 'off', 'TolFun', 1e-6, 'TolX', 1e-6));

param = [];
param.maxIter = 100;      % max number of iterations
param.maxFnCall = 300;  % max number of calling the function
param.relCha = 1e-5;     % tolerance of constraint satisfaction
param.tolPG = 1e-5;      % final objective function accuracy parameter
param.m = 10;

if verbose > 1
    iter_disp = @genericcallback;
else
    iter_disp = [];
end

[sol, loss, ~, ~, msg] = lbfgsb(w0, lb, ub, @local_logistic_loss, [], iter_disp, param); 
if verbose > 0
    fprintf(1, 'Message in svm_logistic_inf: %s\n', msg2str(msg, 'lbfgsb'));
end

sol = sol';

loss = loss / (size(X, 2) * num_task);

function [f, g] = local_logistic_loss(w)    
    F = X' * w;
    f = sum(sum(log(1 + exp(-y .* F))));
        
    M = y ./ (-1 - exp(y .* F));
    g = X * M;    
end

end
