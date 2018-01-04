function X = prox_op_L2ball(A, opts)

%    min_X   0.5*||X - A||_F^2
%     s.t.   ||X_:i||_2 <= opts.radius_ball

  norms = sqrt(sum(A.*A));  
  idx = norms > opts.radius_ball;
  
  X = A;
  if ~isempty(idx) && (length(idx) > 1 || idx > 0)
    X(:,idx) = X(:,idx) ./ repmat(norms(idx)/opts.radius_ball, [size(A,1),1]); 
  end
end
