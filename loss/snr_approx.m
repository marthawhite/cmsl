function snr = snr_approx(Z_recons, Z_noisy) 
% Approximates the Noise-to-signal for pixel valued data WITHOUT needing
% the clean data. This function is mostly necessary for cross validating
% where the clean data is not available. Uses
%    snr = mean pixel value / stddev of pixel value
%
% TODO: Currently not an effective replacement for true snr
%       looking for a way to make this a viable option.
%
% Author: Martha White, University of Alberta, 2012

  if nargin == 0
    snr = 'SNR-Approx';
  elseif nargin == 1 
    fprintf(1, 'Not what you want');
    snr = norm(mean(Z_recons) ./ std(Z_recons));
    % Industry standard: 20 log10 SNR
    %snr = 20 * log10(snr);
    %elseif nargin == 2 % i.e. constrast ratio
    %snr = norm(mean(Z_recons-Z_noisy) ./ std(Z_recons));
  elseif nargin == 2 % i.e. compare to mean and std of clean data
    num_samples = size(Z_recons, 2);
    nsr = sum(abs(mean(Z_recons) - mean(Z_noisy)) + (abs(std(Z_recons) - std(Z_noisy))));
    snr = num_samples/nsr;
  end  
end
