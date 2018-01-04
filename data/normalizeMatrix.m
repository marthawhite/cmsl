Function X = Normalizematrix(X) 
% NORMALIZEMATRIX makes entries in X have variance 1
% Each column of X is a training example.  So divide it by its radius.
  
  num_fea = size(X, 1);
  X = double(X) - repmat(mean(X, 1), num_fea, 1);
  radii = sum(X.^2,1);
  X = X./repmat(sqrt(radii),num_fea,1);
end