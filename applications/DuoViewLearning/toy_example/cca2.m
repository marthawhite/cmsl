function [u,v,cor,err,a,b,phit] = cca2(X,Y)
% uses reduction to pca

[t,n] = size(X);
k = size(Y,2);


[U,S,V] = svd(X,0);
[Q,G,R] = svd(Y,0);
Ztilde = [V*U'; R*Q'];
%SXX = X'*X;
%SYY = Y'*Y;
%SXXih = SXX^(-0.5);
%SYYih = SYY^(-0.5);
%Ztilde = [SXXih*X'; SYYih*Y'];

% pca
%Ktilde = Ztilde*Ztilde';
%[Phi,D] = eig(Ktilde);
[Phi,SS,VV] = svd(Ztilde,0);
s = diag(SS);
%d = diag(D);
d = s.^2;
[dmax,imax] = max(d);
wmax = Phi(:,imax);
a = wmax(1:n);
b = wmax(n+1:end);
phit = wmax'*Ztilde;

% transform
u = V*inv(S)*V'*a*sqrt(2);
v = R*inv(G)*R'*b*sqrt(2);
%u = SXXih*a*sqrt(2);
%v = SYYih*b*sqrt(2);

cor = dmax - 1;
%err = 2*(1 - cor);
err = trace((Ztilde - wmax*phit)'*(Ztilde - wmax*phit));
%norm(X*u)
%norm(Y*v)
