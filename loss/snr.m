function snr = snr(Z, Z_recons) 
% To test (X, X_recons, Y, Y_recons), pass Z-recons = [X_recons; Y_recons]
% Author: Martha White, University of Alberta, 2012

  if nargin < 2
    snr = 'SNR';
  else  
    noise = norm(Z_recons - Z, 'fro')^2;
    signal = norm(Z, 'fro')^2;  
    snr = signal/noise; 
  end
end
