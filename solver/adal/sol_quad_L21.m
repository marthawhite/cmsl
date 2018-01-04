function [X] = sol_quad_L21(A, X0, opts)
% minimize over X (given A and alpha):
%       0.5*|| X - A ||_F^2 + alpha * ||X||_2,1
%
% There is closed form solution!
%
% Input:
%   A: matrix A
%   alpha: weight alpha
%   X0: initial value of X (ignored, just to keep interface unified)
%   opts: options
%
% Output:
%   X: the optimal solution


if nargin ~= 3
    error('sol_quad_L21 requires exactly 3 arguments');
end
  
if opts.reg_wgt < 1e-5
  X = A;  
  return;
end
  
norms = sqrt(sum(A.*A, 2));
idx = find(norms > opts.reg_wgt);
r = zeros(size(norms));
r(idx) = (norms(idx) - opts.reg_wgt) ./ norms(idx);

X = repmat(r, 1, size(A, 2)) .* A;

end
