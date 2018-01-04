function X = prox_op_L21(A, opts)

%    min_X   0.5*||X - A||_F^2 + opts.reg_wgt * ||X||_21

  X = zeros(size(A));
  
  norms = sqrt(sum(A.*A, 2));
  
  reg_wgt = opts.reg_wgt / opts.weight_loss;
  idx = norms > reg_wgt;
  subset = norms(idx);
  
  X(idx,:) = A(idx,:) .* repmat((subset-reg_wgt)./subset, [1, size(A,2)]);

end
