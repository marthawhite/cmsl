function [f,varargout] = barlin(x,A,b,mu)

f = -sum(log(A*x-b))/mu;

if nargout >= 2	% gradient
	varargout{1} = -A'*(1./(A*x-b))/mu;
end
if nargout >= 3	% Hessian
	varargout{2} = A'*diag((A*x-b).^-2)*A/mu;
end
