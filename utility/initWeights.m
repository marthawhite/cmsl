function W = initWeights(sizeW)
% returns an intiial random W, standardized to mean zero
% and variance 1. This is initialized basis according to Lee&Ng's work

    W = randn(sizeW);
    radii = sqrt(sum(W.^2));
    W = W ./ repmat(radii, sizeW(1), 1);       % standardizing

end