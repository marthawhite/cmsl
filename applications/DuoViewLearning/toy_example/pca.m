function V = pca(X)

[V, D] = eig(X'*X);
V = V(:,1);

end
