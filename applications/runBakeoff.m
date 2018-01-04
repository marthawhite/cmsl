function [results, Xinfo, Yinfo, ModelInfo, params, bestParamIndices] = ...
      runBakeoff(loadData, splitData, dataNameMain, AlgNames, opts)
% Returns all the Xtests as an array, unless cross validating, in which case it only returns
% them for the best parameters
%    
% This function tests the objective value, time, and test
% error of the given algorithms over the set of parameters from getParameters
%
% Parameters:
%   loadData: returns the required data for running the experiment, so that it is
%       compatible with the given splitData function. Format of the function is:
%       [X_noisy, Y_noisy] = loadData(dataNameMain, opts) OR
%       [X_noisy, Y_noisy, X_clean, Y_clean] = loadData(dataNameMain, opts).
%   splitData: returns the split for the given index, in case X is returned as
%       a cell array of matrices, corresponding to splits, from loadData.
%       Format of the function is:
%       [X_noisy, Y_noisy] = splitData(X_noisy, Y_noisy, idx, opts)
%
% See below DEFAULTS for optimization options and explanations of those options.
% Note that options are split up into opts.data_params, opts.problem_params and 
% opts.solver_params. There are also other options in top level opts.
%
%   Returns 
%     results.obj_all = zeros(numObjVars,opts.data_params.num_splits,numAlgs); 
%     results.rts_all = zeros(numRuntimeVars,opts.data_params.num_splits,numAlgs); 
%     results.obj = zeros(numObjVars, 2, numAlgs); 
%     results.rts = zeros(numRuntimeVars, 2, numAlgs); 
%     results.err = zeros(numErrVars, 2, numAlgs);
%   where the dimension with 2 includes the means and confidence intervals.
%
% Author: Martha White, University of Alberta, 2012


if nargin < 3
  error('runBakeoff must take at least 3 arguments: functions to loadData, splitData and the data file name.\n');
elseif nargin == 3
  AlgNames =  {'Staged-Alt', 'Alternator', 'Single-View-Col', 'Single-View', 'Separate-View', 'Multi-View'};
end

global col_print;
global sigma_smooth_L11;

MAX_SAMPLES = 2000;  % The maximum number of samples matlab seems to be able to handle with all the algs.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DEFAULT PARAMETERS STARTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DEFAULTS.verbose = 0;
DEFAULTS.output_file = 1;   % default is to output to stdout (i.e. 1)

% If set add_noise function, can include desired options to add_noise
% within data_params; this add noise function is used to aviod having to restore
% huge datasets with different levels of pixel noise (for example)
DEFAULTS.data_params = struct( ...
    'num_splits', 5, ...   % number of times to split the whole datasets into training and testing
    'tl', 100, ...
    'test', MAX_SAMPLES, ...
    'add_noise', [] ...   % additional ability to change noise added (e.g. pixel noise) 
    );

% cv_params.cross_validate 
% == 0 picks one "good" set of parameters, returned by getParameters
% == 1 does cross validation on all params
% == 2 run on all params, pick best result
% If == 2 (i.e. AlgConst.OVER_BEST_PARAMS), then comparison_index specifies
% which error metric to compare them by in result_params.error_fcns
DEFAULTS.cv_params = struct( ...
    'comparison_index', 1, ... 
    'cross_validate', 0, ...   % See AlgConst for cross validate options
    'error_fcn', [], ...       % The function minimized in CV to find best parameters
    'max_samples', 50, ...     % Maximum number of samples to cross validate with
    'num_folds', 5 ...
);

% requires specification of losses, but leaves DEFAULTS value
% for regularizer, L2_wgt, etc., to the algorithms
DEFAULTS.problem_params = struct( ...
    'L1', @smooth_L11, ...  % Default loss is smooth_L11 for X rconstruction
    'L2', @smooth_L11, ...  % Default loss is smooth_L11 for Y reconstruction
     ...%'L2_wgt', 1, ...   % Weight on second loss, to trade-off between reconstruction error for X and Y
    'override_params', 0, ... % if == 1, allows the given values to replace those returned by getParameters     
    'reg_loss', @L21_loss, ...   % Regularization loss
     ...%'reg_wgt', 1e-4,   % Regularization parameter
    'sigma', 0.01 ...     % The smoothness parameter for smooth_L11
    );

% User can modify funTol, maxiter, etc. with solver_params, in case want to
% decrease runtime or increase accuracy. By default, Recover is false (i.e no recovery)
% Possible fields include: 'inner_solver', 'ssl_solver', 'latent_dim', 'num_reps', 'maxiter', 'funTol'
% See individual algorithms for exact parameters (e.g. sparse_coding_convex_primal).
DEFAULTS.solver_params = struct('recover', 0);

% Specify the error functions on the results
% Note that calling an error fcn without parameters should return the name of the error fcn e.g. SNR.
DEFAULTS.result_params = struct( ...
    'error_fcns', {@(X_clean, X_recons, Y_clean, Y_recons)(...
        smooth_L11(X_clean, X_recons) + smoothL11(Y_clean, Y_recons))} ...
    );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DEFAULT PARAMETERS ENDS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 5
  opts = DEFAULTS;
else  
  opts = getOptions(opts, DEFAULTS);
end

sigma_smooth_L11 = opts.problem_params.sigma;

% Get the algorithms to run, with their corresponding parameters
% params is a cell array, that could contain many sets of possible
% parameters for each algorithm is opts.cv_params.cross_validate  == 1
algs = getAlgs(AlgNames);
params = getParameters(AlgNames, (opts.cv_params.cross_validate==AlgConst.OVER_DEFAULT_PARAMS));
if opts.problem_params.override_params
  for a = 1:length(params)
    params{a} = overrideOptions(params{a}, opts.problem_params);
  end  
end
if isempty(algs)
  error('No valid algorithms passed to runBakeoff\n');
end
numAlgs = length(algs);

% Also add any solver or problem parameters to the algorithm parameters
for a = 1:length(params)
  for o = 1:length(params{a})
    params{a}{o} = getOptions(params{a}{o}, opts.problem_params);
    params{a}{o} = getOptions(params{a}{o}, opts.solver_params);
  end 
end

% load the dataset; different experiments will load data different,
% so much provide function for loading experiments
% If function does not provide clean data, then it is simply
% set to empty
X_noisy = [];   Y_noisy = []; X_clean = []; Y_clean = [];
if nargout(loadData) < 4
  [X_noisy, Y_noisy] = loadData(dataNameMain, opts.data_params);
else  
  [X_noisy, Y_noisy, X_clean, Y_clean] = loadData(dataNameMain, opts.data_params);
end
num_ex = size(X_noisy, 2);    % number of examples

% If doing cross validation, take subset of data to pick parameters
if opts.cv_params.cross_validate == AlgConst.CV_PARAMS
  idx = randperm(num_ex);
  numtest = min(opts.cv_params.max_samples, opts.data_params.tl);
  X_cross = X_noisy(:, idx(1:numtest));
  Y_cross = Y_noisy(:, idx(1:numtest));
end

% Variables to record the experimental result
maxNumParams = getMaximumLength(params);
numRuntimeVars = length(IndexConst.RuntimeNames) - double(opts.problem_params.recover == 0);
numErrVars = length(opts.result_params.error_fcns);
numObjVars = length(IndexConst.ObjectiveValueNames);
objs = zeros(opts.data_params.num_splits, numAlgs, maxNumParams, numObjVars);
errs =  zeros(opts.data_params.num_splits, numAlgs, maxNumParams, numErrVars);
runtimes =  zeros(opts.data_params.num_splits, numAlgs, maxNumParams, numRuntimeVars);

% Xinfo contains Xtest_noise, Xtest_clean, Xtest_rec in that order
% See IndexConst for all the indexing constants
Xinfo = cell(opts.data_params.num_splits, 1);
Yinfo = cell(opts.data_params.num_splits, 1);
ModelInfo = cell(opts.data_params.num_splits, 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RUN EXPERIMENT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for idx_split = 1 : opts.data_params.num_splits   % for each split of training/test set
  
  % Let calling experiment chose how to split data (done differently for duo-view, semisupervised)
  % And different for classification and regression
  [Xi, Yi, Xi_clean, Yi_clean, Xtest_noise, Ytest_noise, Xtest_clean, Ytest_clean] = ...
      splitData(X_noisy, Y_noisy, X_clean, Y_clean, idx_split, opts);
  
  Xinfo{idx_split}{IndexConst.CLEAN_IND} = Xtest_clean;
  Xinfo{idx_split}{IndexConst.NOISE_IND} = Xtest_noise;
  Yinfo{idx_split}{IndexConst.CLEAN_IND} = Ytest_clean;
  Yinfo{idx_split}{IndexConst.NOISE_IND} = Ytest_noise;
  
  fprintf(1,'*********************************************');
  fprintf(1, '*********************************************************\n');
  fprintf(1,'Running split = %d for tl = %d and n = %d\n', idx_split, opts.data_params.tl, size(Xi,1));
  fprintf(1,'*********************************************');
  fprintf(1, '*********************************************************\n');
    
  Xinfo{idx_split}{IndexConst.RECON_IND} = cell(numAlgs); 
  Yinfo{idx_split}{IndexConst.RECON_IND} = cell(numAlgs);   
  for a = 1:numAlgs
    alg = algs{a};
    % Only cross validate once on subset of data
    if idx_split == 1 && opts.cv_params.cross_validate == AlgConst.CV_PARAMS
      cv_start = cputime;   
      [minIndex, minParams] = crossValidate(params{a}, alg, X_cross, Y_cross, opts.cv_params);
      params{a} = {minParams};
      if opts.verbose >= VerboseConst.MAIN_SCRIPT
        fprintf(1, 'runBakeoff -> Took %d seconds to find cv parameters for alg = %s for datatset = %\n', ...
                cputime-cv_start, AlgNames{a}, dataNameMain);
      end    
    end    
    Xinfo{idx_split}{IndexConst.RECON_IND}{a} = cell(length(params{a})); 
    Yinfo{idx_split}{IndexConst.RECON_IND}{a} = cell(length(params{a})); 
    Xtest_rec = []; Ytest_rec = []; Ball = [];Phiall = []; Wall = [];
    for o = 1:length(params{a})
      algopts = params{a}{o};
      
      if opts.verbose >= VerboseConst.MAIN_SCRIPT_DETAILED
        fprintf(1,'***************************************************************\n');
        fprintf(1,'Running split = %d, algorithm = %s, \n', idx_split, AlgNames{a});
        fprintf(1, 'algopts = %s\n', struct2str(algopts));
        fprintf(1,'***************************************************************\n');
      end      
      algopts.verbose = opts.verbose;
      
      % Call the current alg with the current params; if no recovery
      % test error reported transductively
      if opts.problem_params.recover
        [X_recons Y_recons pobj runtime B W Phi] = alg(Xi, Yi, algopts);
        % Re-adjust reconstruction      
        [X_recons Y_recons Xtest_recons Ytest_recons B W Phi pobj] = ...
            adjust_recovery(Xi, Yi, Xtest_noise, Ytest_noise, B, W, Phi, opts.solver_params);        
    		ModelInfo{idx_split}{IndexConst.B_IND}{a}{o} = B;   
    		ModelInfo{idx_split}{IndexConst.W_IND}{a}{o} = W;   
    		ModelInfo{idx_split}{IndexConst.PHI_IND}{a}{o} = Phi;        
        objs(idx_split, a, o, IndexConst.REG_IND) = opts.problem_params.reg_loss(Phi);        
      else
        [X_recons Y_recons pobj runtime] = alg(Xi, Yi, algopts); 
        Xtest_recons = X_recons;
        Ytest_recons = Y_recons;
        runtime = runtime(1);
      end
      % Save all solutions for that algorithm
      Xinfo{idx_split}{IndexConst.RECON_IND}{a}{o} = Xtest_recons;
      Yinfo{idx_split}{IndexConst.RECON_IND}{a}{o} = Ytest_recons;      
      
      % Record the objective values
      objs(idx_split, a, o, IndexConst.OBJ_IND) = pobj;
      objs(idx_split, a, o, IndexConst.L1_IND) = opts.problem_params.L1(Xi, X_recons);
      objs(idx_split, a, o, IndexConst.L2_IND) = opts.problem_params.L2(Yi, Y_recons);
      
      % Record the error values
      for err_f = 1:numErrVars
        errs(idx_split, a, o, err_f) = opts.result_params.error_fcns{err_f}(...
                            Xtest_clean, Xtest_recons, Ytest_clean, Ytest_recons); 
      end

      runtimes(idx_split, a, o, :) = runtime;
      if opts.verbose >= VerboseConst.MAIN_SCRIPT_DETAILED
        col_print(['\nobj = %g, L1-xrecon = %g, L2-yrecon = %g, clean-xrecon = %g, '... 
                   'clean-yrecon = %g, time = %g, error1 = %g\n\n'], ...
                  pobj, objs(idx_split, a, o, IndexConst.L1_IND), objs(idx_split, a, o, IndexConst.L2_IND), ...
                  opts.problem_params.L1(Xi_clean,X_recons), opts.problem_params.L2(Yi_clean,Y_recons),...
                  runtime(IndexConst.TRAINING_TIME_IND), errs(idx_split, a, o, 1)); 
      end
    end     % parameters
  end     % algorithm
end     % data split index


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% OUTPUT RESULTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf(1, ['Warning! The objectives are sometimes not comparable because the ' ...
            'algorithms are optimizing different functions or have different '...
            'parameters. They are mostly for sanity checking!\n']);

% Open output file, appending to previous file 
if isempty(opts.output_file) || isnumeric(opts.output_file)
  fileId = 1;
else 
  fileId = fopen(opts.output_file, 'a+');
  if (fileId == -1)
    fileId = 1;
  end    
end
fprintf(fileId, '**********************************************\n');
fprintf(fileId, 'On date = %s, num_splits = %d, tl = %d on datatset = %s\n\n', ...
        date, opts.data_params.num_splits, opts.data_params.tl, dataNameMain);

% Output results
sqrtN = sqrt(opts.data_params.num_splits);
results = [];
separator = '    ';
algStringLength = getMaximumLength(AlgNames) + 3;
columnLen = length('6.544e+01 +/- 7.305e+00');

fprintf(fileId, '%-*s', algStringLength, 'Algorithm');
for err_f = 1:numErrVars 
  funcName =  opts.result_params.error_fcns{err_f}();
  fprintf(fileId, [separator '%-*s'], columnLen, funcName);
end
for l = 1:numObjVars
  if opts.problem_params.recover == 0 && l == IndexConst.REG_IND, continue; end
  fprintf(fileId, [separator '%-*s'], columnLen, IndexConst.ObjectiveValueNames{l});
end
fprintf(fileId, [separator '%-*s'], columnLen, 'Runtime Training');
if opts.problem_params.recover ~= 0
  fprintf(fileId, [separator '%-*s'], columnLen, 'Runtime Recovery');
end
fprintf(fileId, '\n\n');

results.obj_all = zeros(numObjVars,opts.data_params.num_splits,numAlgs); 
results.rts_all = zeros(numRuntimeVars,opts.data_params.num_splits,numAlgs); 
results.obj = zeros(numObjVars, 2,numAlgs); 
results.rts = zeros(numRuntimeVars, 2,numAlgs); 
results.err = zeros(numErrVars,2,numAlgs);
numParams = max(cellfun('length', params));
results_params = [];
if numParams > 1 && (opts.cv_params.cross_validate ~= AlgConst.OVER_BEST_PARAMS)
	results_params.obj_all = zeros(numObjVars,opts.data_params.num_splits,numAlgs,numParams); 
	results_params.rts_all = zeros(numRuntimeVars,opts.data_params.num_splits,numAlgs,numParams); 
	results_params.obj = zeros(numObjVars, 2,numAlgs,numParams); 
	results_params.rts = zeros(numRuntimeVars, 2,numAlgs,numParams); 
	results_params.err = zeros(numErrVars,2,numAlgs,numParams);  
end

% Assume every alg has first parameters as best;
% OVER_BEST_PARAMS is selected, then these indices will
% be changed; for CV_PARAMS and OVER_DEFAULT, length(params) = 1
bestParamIndices = ones(1,numAlgs);
for a = 1:numAlgs
  best_mean_err = inf; % Only used if opts.cv_params.cross_validate == AlgConst.OVER_BEST_PARAMS
  numParams = length(params{a});
  for o = 1:numParams
    mean_obj = mean(objs(:, a, o, :),1);
    mean_rts = mean(runtimes(:, a, o,:),1);
    mean_err = mean(errs(:,a, o,:),1);
    std_obj = std(objs(:, a, o,:),0,1)/sqrtN;
    std_rts = std(runtimes(:, a, o,:),0, 1)/sqrtN;
    std_err = std(errs(:, a, o,:),0,1)/sqrtN;
    
    if numParams <= 1
      results.obj_all(:,:,a) = squeeze(objs(:, a, o, :))'; 
      results.rts_all(:,:,a) = squeeze(runtimes(:, a, o, :))'; 
      results.obj(:,:,a) = [squeeze(mean_obj) squeeze(std_obj)]; 
      results.rts(:,:,a) = [squeeze(mean_rts) squeeze(std_rts)];
      results.err(:,:,a) = [squeeze(mean_err) squeeze(std_err)];
    elseif opts.cv_params.cross_validate == AlgConst.OVER_BEST_PARAMS
      if mean_err(opts.cv_params.comparison_index) < best_mean_err
        results.obj_all(:,:,a) = squeeze(objs(:, a, o, :))'; 
        results.rts_all(:,:,a) = squeeze(runtimes(:, a, o, :))'; 
        results.obj(:,:,a) = [squeeze(mean_obj) squeeze(std_obj)]; 
        results.rts(:,:,a) = [squeeze(mean_rts) squeeze(std_rts)];
        results.err(:,:,a) = [squeeze(mean_err) squeeze(std_err)];
        best_mean_err = mean_err(opts.cv_params.comparison_index);
        bestParamIndices(a) = o;
      end
      if o < numParams, continue; end
    else 
      results.obj_all(:,:,a,o) = squeeze(objs(:, a, o, :))'; 
      results.rts_all(:,:,a,o) = squeeze(runtimes(:, a, o, :))'; 
      results.obj(:,:,a,o) = [squeeze(mean_obj) squeeze(std_obj)]; 
      results.rts(:,:,a,o) = [squeeze(mean_rts) squeeze(std_rts)];
      results.err(:,:,a,o) = [squeeze(mean_err) squeeze(std_err)];
    end
  
    vars = [];
    for num = 1:numErrVars
      vars = [vars results.err(num,:,a)];
    end  
    for num = 1:numObjVars
      if opts.problem_params.recover == 0 && num == IndexConst.REG_IND, continue; end      
      vars = [vars results.obj(num,:,a)];
    end      
    vars = [vars results.rts(IndexConst.TRAINING_TIME_IND,:,a)];
    if opts.problem_params.recover == 1
      vars = [vars results.rts(IndexConst.RECOVERY_TIME_IND,:,a)];
    end    
    
    fprintf(fileId, '%-*s', algStringLength, AlgNames{a});
    ind = 1;
    while ind < length(vars)
      fprintf(fileId, [separator '%.3i +/- %.3i'],vars(ind), vars(ind+1));
      ind = ind+2;
    end    
    fprintf(fileId, '\n');
  end
  if ~isempty(results_params)
    results_params.obj_all(:,:,a,o) = results.obj_all(:,:,a); 
    results_params.rts_all(:,:,a,o) =  results.rts_all(:,:,a)
    results_params.obj(:,:,a,o) =  results.obj(:,:,a) ; 
    results_params.rts(:,:,a,o) = results.rts(:,:,a);
    results_params.err(:,:,a,o) = results.err(:,:,a);
  end
end
if ~isempty(results_params)
  results = results_params;
end

for a = 1:numAlgs
  cprintf(VerboseConst.MAIN_SCRIPT, '\nBest parameters for %s are %s\n',...
          AlgNames{a}, struct2str(params{a}{bestParamIndices(a)}));  
end

if (fileId ~= 1)
  fclose(fileId);
end

end
