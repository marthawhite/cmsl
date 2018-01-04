%TEST
clear;clc;

initialize_testing_paths;

opts_array =  {struct('reg_wgt', 0.01, 'gamma', 0.1), struct('reg_wgt', 0.1, 'gamma', 0.1), ...
	     struct('reg_wgt', 0.01, 'gamma', 1), struct('reg_wgt', 0.1, 'gamma', 1)};

% First test if can add new options to all options
added_opts =  struct('maxiter', 500, 'num_reps', 5);
augmented_array = addToAllOptions(opts_array, added_opts);
fieldname_array = {'reg_wgt','gamma','maxiter','num_reps'};
if length(augmented_array) ~= length(opts_array) || ...
      ~all(isfield(augmented_array{1}, fieldname_array))
  fprintf(1, 'addToAllOptions failed.\n');
else
  fprintf(1, 'addToAllOptions succeeded.\n');
end  

% Second test if can overrideOptions
reg_wgts = [0.5 1 10];
num_gamma = 2;
override_opts =  struct('reg_wgt', reg_wgts);
override_array = overrideOptions(opts_array, override_opts);
fieldname_array = {'reg_wgt','gamma'};
if length(override_array) ~= num_gamma*length(reg_wgts)
  fprintf(1, ['overrideOptions failed, number of options %d in opts_array'...
      'incorrect, not equal to %d.\n'], length(override_array), num_gamma*length(reg_wgts));
  for i = 1:length(override_array)
    override_array{i}
  end  
elseif ~all(isfield(override_array{1}, fieldname_array))
  fprintf(1, 'overrideOptions failed, removed a field incorrectly from opts_array.\n');
elseif override_array{1}.reg_wgt ~= 0.5
  fprintf(1, 'overrideOptions failed, did not correctly override reg_wgt.\n');
else
    fprintf(1, 'overrideOptions succeeded.\n');
end 


