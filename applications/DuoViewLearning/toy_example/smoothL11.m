function f = smoothL11(x, sigma)

		idx = abs(x) < sigma;
    Z = abs(x) - sigma/2;
    Z(idx) = x(idx).^2 * (0.5/sigma);
    
    f = sum(sum(Z)); 

end