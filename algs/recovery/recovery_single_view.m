function [B, W, Phi, runtime] = recovery_single_view(X_recons, Y_recons)

  start_time = cputime;
  % Now recover B, W and Phi from Zhat
  Zhat = [X_recons; Y_recons];
  [Usvd,Sigma,V] = svd(Zhat, 'econ'); 
  U = Usvd;
  Phi = Sigma * V';
  B = U(1:size(X_recons,1), :);
  W = U(1+size(X_recons,1):end, :);
  runtime = cputime - start_time;  

end
