  function [f, g] = pobj_adal(X, L1, L2, A, reg_wgt)
    X = reshape(X, size(A));
    [f1, g1] = L1(X);
    [f2, g2] = L2(X);
    f = f1 + reg_wgt * f2;
    g = g1 + reg_wgt * g2;
    g = g(:);
  end