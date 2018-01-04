function [x,P,flag] = feastarteq(Aeq,beq,x0)
% project x0 onto Aeq*x = beq

TOL = 1e-8;

if rank(full(Aeq)) < size(Aeq,1) 
	warning('equality constraint matrix not full rank'); 
	flag = 1;
	x = x0;
	P = NaN;
	return
end

flag = 0;
P = (Aeq*Aeq') \ Aeq;
%Proj = eye(length(x0)) - Aeq'*P;	% projection matrix
x = x0;
res = Aeq*x - beq;
if max(abs(res)) > TOL
	y0 = (Aeq*Aeq') \ res;
	d0 = -Aeq'*y0;
	x = x0 + d0;
end

