function [W, loss] = svm(X,y,alpha)
% X is n x t with t the number of samples
% y is 1 x t in {-1,1}
% if alpha = 1/(2beta) is provided, then the
% soft margin SVM is run
% Notation from Xu and Schuurmans paper

X = X';
y = y';
[t,n] = size(X);
if (nargin < 3)
	H = [eye(n,n) zeros(n,1); zeros(1,n+1)];
	f = zeros(n+1,1);
	
	% quadprog
	A = -[repmat(y,1,n).*X y];
	b = -ones(t,1);
    warning('off', 'all');
	[Wb,qploss,flag] = quadprog(H,f,A,b,[],[],[],[],[],optimset('Display', 'off'));
    
    loss = qploss;
	W = Wb(1:n);
else    
    K = X*X';
	H = K.*(y*y')/alpha;
	f = -ones(t,1);
	
	% quadprog
	lb = zeros(t,1);
	ub = ones(t,1);
	[lambda,qploss,flag] = quadprog(H,f,[],[],[],[],lb,ub, [],optimset('Display', 'off'));
	
    loss = -qploss;
	W = X'*(lambda.*y)/alpha;    
end
W = W';

