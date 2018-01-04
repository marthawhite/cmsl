function [R] = invhessmult(V,Y,S,Rho,H0,inds,m)
% implicit multiplication of V by limited memory inverse Hessian approximation

[t,n] = size(V);
Alpha = zeros(length(inds),n);
gamma = 1;

Q = V;
for j = length(inds):-1:1
	i = inds(j);
	Alpha(i,:) = Rho(i)*S(:,i)'*Q;
	Q = Q - Y(:,i)*Alpha(i,:);
end
if length(inds) == m
	gamma = S(:,inds(1))'*Y(:,inds(1)) / (Y(:,inds(1))'*Y(:,inds(1)));
end
R = gamma*H0*Q;
for j = 1:length(inds)
	i = inds(j);
	Beta = Rho(i)*Y(:,i)'*R;
	R = R + S(:,i)*(Alpha(i,:)-Beta);
end
