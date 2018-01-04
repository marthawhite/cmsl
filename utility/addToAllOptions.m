function opts_array = addToAllOptions(opts_array, added_opts) 
% opts_array is an array of options, where this function
% adds all the fields in added_opts to each opts in the array
% If any opts in opts_array already contains the field, it does not
% override the field (use overrideOptions for that functionality).
%
% e.g. s = addToAllOptions(...
%  {struct('reg_wgt', 0.01, 'gamma', 0.1), struct('reg_wgt', 0.1, 'gamma', 0.1),…
%	     struct('reg_wgt', 0.01, 'gamma', 1), struct('reg_wgt', 0.1, 'gamma', 1)},…
%	 struct('maxiter', 500));
%
% returns:
% s =  {struct('reg_wgt', 0.01, 'gamma', 0.1, 'maxiter', 500), ...
%       struct('reg_wgt', 0.1, 'gamma', 0.1 'maxiter', 500),…
%       struct('reg_wgt', 0.01, 'gamma', 1, 'maxiter', 500),...
%       struct('reg_wgt', 0.1, 'gamma', 1, 'maxiter', 500)},
% 
% If opts_array is a struct full of structs, then each struct has
% added_opts added susing getOptions
%
% Author: Martha White, University of Alberta, 2012
  
  if isempty(opts_array) || isempty(added_opts)
	  return;
  end
  
  if isstruct(opts_array)
    field_names = fieldnames(opts_array)';    
    for i=1:length(field_names)
     f = field_names{i};
     fld = getfield(opts_array, f);
     if isstruct(fld)
       opts_array = setfield(opts_array, f, getOptions(fld, added_opts));
     end
    end  
    return;
  end
  
  for i = 1:length(opts_array)
    opts_array{i} = getOptions(opts_array{i}, added_opts);
  end

end
