function opts = overrideOptions(opts, opts2) 
% Enables override of variables within opts, if opts2 also has them 
% If opts is a cell array of many options, replace for each
  
  if nargin < 2
		return;
  elseif ~iscell(opts)
    opts = replaceFields(opts, opts2);
    return;
  end
  
  % If opts2 contains a list where opts contains only a scalar,
  % then add a bunch of duplicates to params for that repeat
  opts2_fields = fieldnames(opts2);
  for i = 1:length(opts2_fields)
    fname = opts2_fields{i};
    f_opts2 = opts2.(fname);
    is_list = 0;
    for i = 1:length(opts)
      if isfield(opts{i}, fname)
        f_opts1 = opts{i}.(fname);
        if strcmp(class(f_opts1),class(f_opts2)) && length(f_opts2) > 1
          %  if isfield(opts2, 'reg_wgt') && length(opts2.reg_wgt) > 1
          is_list = 1;
          break;
        end
      end
    end
    if is_list
      newParams = [];
      for i = 1:length(opts)
        for a = 1:length(f_opts2)
          opts{i}.(fname) = f_opts2(a);
          newParams = [newParams {opts{i}}];
        end  
      end  
      opts2 = rmfield(opts2, fname);
      opts = newParams;
    end
  end
  
  % Remove any params that become redudant from replacement  
  newParams = [];
  for i = 1:length(opts)
    opts{i} = replaceFields(opts{i}, opts2);
    if (i < 2 || alreadyContainsParams(newParams, opts{i}, i-1) == 0)
      newParams = [newParams {opts{i}}];
    end
  end      
  opts = newParams;
  
  function opts = replaceFields(opts, opts2) 
    fields = fieldnames(opts);
    for i2=1:length(fields),
      f = fields{i2};
      if (isfield(opts2,f))
        opts.(f) = opts2.(f);
      end
    end
  end
  
  function bool = alreadyContainsParams(paramsArray, params, finalIndex)
    for i2 = 1:min(finalIndex, length(paramsArray))
      if (isequal(params,paramsArray{i2}))
        bool = 1;
        return;
      end
    end
    bool = 0;
  end    
  
end
