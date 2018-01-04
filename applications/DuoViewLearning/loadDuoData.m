function [X_noisy, Y_noisy, X_clean, Y_clean] = loadDuoData(dataNameMain, opts)

X = []; Y = []; X1 = []; X2 = []; Y1 = []; Y2 = [];
if isfield(opts, 'nfeatures') && opts.nfeatures > 0
  [X Y X1 X2 Y1 Y2] = generate_data(opts.nfeatures, opts.nfeatures, ...
                                    opts.tl+opts.test, opts);
else 
  load(dataNameMain);  % This will fill X and Y
end 

X_noisy = X; Y_noisy = Y; X_clean = X1; Y_clean = Y1;

% Add any extra noise that applications might want (e.g. pixel noise)
if ~isempty(opts.add_noise)
  [X_noisy, Y_noisy] = opts.add_noise(X_noisy, X_clean, Y_noisy, Y_clean, opts);
end

end
