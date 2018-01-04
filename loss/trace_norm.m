function [val, G] = trace_norm(A) 
% Computes the trace norm on mxn A
% 
%	||A||_tr = sum_i^{min(m,n)} sigma_i
%
% Author: Martha White, University of Alberta, 2011

  %if all(abs(A) < 1e-2)
  %  val = 0;
  %  G = zeros(size(A));
  %  return;
  %end  
  if nargout < 2
    s = svd(A, 'econ');
    val = sum(s);
  else  
    [U, S, V] = svd(A, 'econ');
    val = sum(S(:));  % diag allows S to be rectangular
    G = U*V'; % subgradient
  end

end
