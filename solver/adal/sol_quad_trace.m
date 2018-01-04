function [X] = sol_quad_trace(A, X0, opts)
% minimize over X (given A and reg_wgt):
%       0.5*|| X - A ||_F^2 + reg_wgt * ||X||_tr
%
% There is a closed form solution!
%
% Input:
%   A: matrix A
%   X0: initial value of X (ignored, just to keep interface unified)
%   opts: options
%
% Output:
%   X: the optimal solution
%
% Author: Yaoliang Yu, University of Alberta, 2012


if nargin < 3
    error('sol_quad_trace requires at least 3 arguments: A, X0 and opts. X0 is not used.');
end

[U S V] = svd(A, 'econ');
s = max(0, diag(S) - opts.reg_wgt);
X = (U .* repmat(s', size(U, 1), 1)) * V';

end
