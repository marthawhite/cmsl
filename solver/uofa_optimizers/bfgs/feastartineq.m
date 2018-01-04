function [x,flag] = feastartineq(A,b,x0)
% checks if x0 feasible: A*x >= b
% if x0 infeasible, computes feasible start with LP

TOL = 1e-14;
flag = BFGS_ERRORS.SUCCESS;
x = x0;
[k,n] = size(A);

% check if already feasible
if all(A*x >= b), return; end

% if x0 infeasible, compute feasible start
ff = [zeros(n,1); -1];
AA = [-A ones(k,1); zeros(1,n) 1];
bb = [-b; 1/n];

% linprog
[z,val,lflag] = linprog(ff,AA,bb,[],[],[],[],[],optimset('Display','off'));

% cplex
%ctype = repmat('L',size(AA,1),1);
%vtype = repmat('C',size(AA,2),1);
%[z,val,lflag] = cplexmex(1,[],ff,AA,bb,ctype,[],[],vtype);

flag = lflag - 1;
if isempty(z)
    warning('feastartineq -> linprog could not find a feasible start!');
    x = x0;
else    
    x = z(1:end-1);
end

