function [C,H] = pca(X, latent_dim)


[U,S,V] = svd(X,0);
s = diag(S);
umax = U(:,1:latent_dim);
smax = s(1);
vmax = V(:,1:latent_dim);
C = umax;
H = smax*vmax';
err = sum(s(2:end).^2);
