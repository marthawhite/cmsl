function filename = getDataNameDuo(n, c, t, opts)

  global data_dir;

  DEFAULTS = getDefaultDataParams();

  if nargin < 3
    opts = DEFAULTS;
  else  
    opts = getOptions(opts, DEFAULTS);
  end

  if opts.shifted == 0
    filename = sprintf([data_dir '/DuoView/duoRPCA_n=%d_c=%d_t=%d_rk=%d_rto=%d_mag=%d.mat'], ...
                       n, c, t, opts.rk, opts.ratio, opts.mag);
  else
    filename = sprintf([data_dir '/DuoView/duoRPCA_n=%d_c=%d_t=%d_rk=%d_ratio=%d_mag=%d_sft=1.mat'], ...
                       n, c, t, opts.rk, opts.ratio, opts.mag); 
  end
end
