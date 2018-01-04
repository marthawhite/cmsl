function [minIndex,minParams] = crossValidate(params, alg, X, Yl,opts)
% CROSSVALIDATE finds the best parameters for the given fcn
% based on given data X and Yl, where in supervised or semisupervisd case
% X = [Xl; Xu] and Yl labeled targets and in multi-view case,
% X and Yl are the two views.
%
% The CV is done transductively, since many competing algorithms
% can only be transductive. TODO: extend CV to using models
%
% Author: Martha White, University of Alberta, 2012

% DEFAULTS.error_fcn: Function for comparing optimiality of solution
% error_fcn(X, X_recons, Y, Y_recons) for multiview
% and error_fcn(Y_test, Y_pred) for (semi-)supervised  
DEFAULTS.error_fcn = @(X, X_recons, Y, Y_recons)(-snr_approx([X_recons; Y_recons]));   
DEFAULTS.num_folds = 10;   % Number of folds
DEFAULTS.verbose = 0;   % Whether to print out progress information
DEFAULTS.timeout = 300;   % No algorithm can take more than timeout seconds for each parameter

if nargin < 5
    opts = DEFAULTS;
else
    opts = getOptions(opts, DEFAULTS);
end

num_params = length(params);
if num_params == 1
    minIndex = 1;
    minParams = params{1};
    return;
end
   
% Save original params as they get slightly modified when running looo
original_params = params;

% Determine if we are doing (semi)supervised or multiview learning
tl = size(Yl,2);
err_all = zeros(num_params,1);
is_multiview = 1;
Xl = X(:, 1:tl);
Xu = [];
if nargin(opts.error_fcn) < 4
  is_multiview = 0;
  Xu = X(:, (tl+1):end);
elseif tl ~= size(X, 2)
    error('crossValidate -> error function for multiview learning, but two views have different amount of samples\n');
end

% Create the folds; if tl < num_folds, cap num_folds
p = randperm(tl); 
opts.num_folds = min(opts.num_folds,tl);
numInEachFold = floor(tl/opts.num_folds);
folds = [];
for i = 1:opts.num_folds
    folds = [folds {p((i-1)*numInEachFold+1:i*numInEachFold)}];
end    

if opts.verbose >= VerboseConst.BASIC_CV, fprintf(1, '\nStarting cross validation...\n'); end

% Cross validate with leave-one-out
for s=1:opts.num_folds
  testInd = folds{s};
  trainInd = [];
  for i = 1:opts.num_folds
      if i == s
          continue;
      end
      trainInd = [trainInd folds{i}];
  end    
  Y_t = Yl(:,testInd);
  Yl_c = Yl(:,trainInd);
  Xl_c = Xl(:,trainInd);
  Xu_c = [Xl(:,testInd) Xu];
  if is_multiview
    % Currently ignores test data; will change when add CV with models
  	for i=1:num_params
      params{i}.verbose = opts.verbose;
      params{i}.timeout = opts.timeout;
  	  [X_learned Y_learned] = alg(Xl_c, Yl_c, params{i});
      err = opts.error_fcn(Xl_c, X_learned,Yl_c, Y_learned);
  	  err_all(i) = err_all(i) + err;
      if opts.verbose >= VerboseConst.DETAILED_CV, 
        fprintf(1, 'fold = %u, error = %g, opts = %s\n', s, err, struct2str(params{i})); 
      end      
   end
  else
    offset = length(trainInd)+1;
    for i=1:num_params
      params{i}.verbose = opts.verbose;
      params{i}.timeout = opts.timeout;
  	  [X_learned Y_learned]= alg([Xl_c Xu_c], Yl_c, params{i});
      err = opts.error_fcn(Y_t,Y_learned(:,offset:offset+length(testInd)));
  	  err_all(i) = err_all(i) + err;
      if opts.verbose >= VerboseConst.DETAILED_CV, 
        fprintf(1, '\n****** fold = %u, error = %g, opts = %s\n\n', s, err, struct2str(params{i})); 
      end
  	end
  end
  if opts.verbose >= VerboseConst.BASIC_CV && opts.verbose < VerboseConst.DETAILED_CV 
    fprintf(1, '*'); 
  end
end
if opts.verbose >= VerboseConst.BASIC_CV, fprintf(1, '\n'); end

[minErr,minIndex] = min(err_all);
minIndex = minIndex(1);
minParams = original_params{minIndex};

end