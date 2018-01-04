function [X_noisy,Y_noisy] = add_pixel_noise(X_noisy, X_clean, Y_noisy, Y_clean, opts)

% The percentage noise for X and Y;
% If percentage noise = 0, then sets to clean data
% If percentage noise = -1, then does not change the noise
DEFAULTS.percentage_ones = [0 0];

if nargin < 5
  opts = DEFAULTS;
else
  opts = getOptions(opts, DEFAULTS);
end

X_noisy = add_noise(X_noisy, X_clean, opts.percentage_ones(1));
Y_noisy = add_noise(Y_noisy, Y_clean, opts.percentage_ones(2));

  function Z_noisy = add_noise(Z_noisy, Z_clean, percentage)
    if percentage >= 0 && percentage < 1
      Z_noisy = Z_clean;
      if percentage > 0
        idx = randperm(numel(Z_noisy));
        messed_up = floor(numel(Z_noisy)*percentage);
        Z_noisy(idx(1:messed_up)) = 1;  
      end
    end
  end

end

