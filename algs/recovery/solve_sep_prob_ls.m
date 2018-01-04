% Solve max_b  b' A b
% s.t. b' Sigma b = 0 and b'b = 1

function [basis, t2 t3 t4] = solve_sep_prob_ls(A, Sigma, rho_init, param)

  % First line search over rho  
   tt = cputime;
%   max_rho = 100;
%   [rho, f, ~, ~] = fmin_BFGS(@eig_search, rho_init, ...
%                              [1 -1], [-max_rho -max_rho], [], [], 1e-5, 200);

%   [rho, f, iter, num_call, msg] = lbfgsb(rho_init, -max_rho, max_rho, ...
%                                           @eig_search, [], [], param);
                           
  [rho, f, iter, feval, flag] = pbm(@eig_search, rho_init, -inf, inf, param);                         

%   fprintf('pbm in solve_sep_prob_ls took feval: %d\n', feval);

  t2 = cputime - tt;  warning on;  
  
  % Then recover the optimal b
  tt = cputime;                                            
  G = null_space(diag(rho*Sigma+f) - A);
  t3 = cputime - tt;
  
  M = G'* (repmat(Sigma, 1, size(G, 2)).*G);
  tt = cputime;
  [V m] = eigs(M, 1, 'SM');   % get a vector V in the null spae of M
                              % using eigs for numerical robustness (M >= 0)
  t4 = cputime - tt;
  basis = G * V;
  
%   if abs(m) > 1e-2
%     cprintf('r', 'warning in solve_sep_prob_ls: %g\n', m);
%   end
  
  function [f, g] = eig_search(rho)
    [V f] = eigs(A - rho *diag(Sigma), 1, 'LA'); % greatest algebraic eigenvalue
%     [V,f] = laneig(Q - rho *diag(Sigma), 1, 'AL'); % slow here
    g = -(Sigma'*(V.^2));
  end
end
