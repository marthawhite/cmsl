function opts = getDuoViewParams(percentage_ones)
% If percentage_ones is not given, does not add pixel noise
% If pixel noise added, then sets SNR error function based on
% percentage_ones

opts = [];
  
% Negative of approximate SNR for cross validation; approximate SNR, because
% do not have clean data for cross validation.
opts.cv_params.error_fcn =  @(X, X_recons, Y, Y_recons)(-snr_approx([X_recons; Y_recons],[X;Y]));

opts.data_params = [];
if nargin >= 1
    opts.data_params.add_noise = @add_pixel_noise;
    opts.data_params.percentage_ones = percentage_ones;
end

opts.problem_params.recover = false;

% Only report SNR on Y if all noise removed from X
opts.result_params = [];
if ~isfield(opts.data_params, 'percentage_ones') || opts.data_params.percentage_ones(1) ~= 0
  opts.result_params.error_fcns = {@duoViewSNR1};
else
  opts.result_params.error_fcns = {@duoViewSNR2};
end
opts.result_params.error_fcns = [opts.result_params.error_fcns {@duoViewClean1, @duoViewClean2}];

function s = duoViewClean1(X_clean, X_recons, Y_clean, Y_recons)
  if nargin < 2
    s = 'Clean X Error';
  else
    s = norm(X_clean - X_recons, 'fro')/size(X_clean,2);
  end  
end

function s = duoViewClean2(X_clean, X_recons, Y_clean, Y_recons)
  if nargin < 2
    s = 'Clean Y Error';
  else
    s = norm(Y_clean - Y_recons, 'fro')/size(Y_clean,2);
  end  
end

function s = duoViewError(X_clean, X_recons, Y_clean, Y_recons)
  if nargin < 2
    s = 'L2 Error';
  else
    s =  norm([X_clean; Y_clean] - [X_recons; Y_recons], 'fro')/size(X_clean,2);
  end  
end

function s = duoViewSNR1(X_clean, X_recons, Y_clean, Y_recons)
  if nargin < 2
    s = '-SNR';
  else
    s = -snr([X_clean;Y_clean], [X_recons;Y_recons]);
  end  
end

function s = duoViewSNR2(X_clean, X_recons, Y_clean, Y_recons)
  if nargin < 2
    s = '-SNR (on Second View)';
  else
    s = -snr(Y_clean, Y_recons);
  end  
end

end
