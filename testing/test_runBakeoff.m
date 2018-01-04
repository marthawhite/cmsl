% TEST

% Tests runBakeoff with different result options
clear;clc;

global proj_dir;
initialize_testing_paths;
addpath([proj_dir '/applications/DuoViewLearning']);

% Test with fast algorithms
AlgNames = { 'Multi-View', 'Single-View'};
numAlgs = length(AlgNames);

% Use noiseless data to make faster and to test
% the add noise aspect of runBakeoff
opts = getDuoViewParams([0 0]);
opts.verbose = VerboseConst.MAIN_SCRIPT;

% Make data small
gopts = [];
gopts.rk = 5; gopts.ratio = 5; gopts.mag = 10; 
t = 5000; % Total points in dataset
num_features = 20;
fileName = getDataNameDuo(num_features,num_features,t,gopts);
if exist(fileName) == 0
  generateDataDuo(num_features,num_features,t,gopts);
else
  fprintf(1, 'Loading existing file %s\n', fileName);
end 
opts.data_params.test = 0;
opts.data_params.tl = 200;
opts.data_params.num_splits = 5;

% Add a timeout to ensure solvers do not take too long
opts.solver_params.timeout = 100;
opts.solver_params.recover = 0;

% Ensure that can run parameters over default params
fprintf(1, '++++++++++++++++++++++++++++++++ Running test for OVER_DEFAULT_PARAMS...\n\n');
opts.cv_params.cross_validate = AlgConst.OVER_DEFAULT_PARAMS;
[results, Xinfo, Yinfo, ModelInfo, params, bestParamIndices] = ...
    runBakeoff(@loadDuoData, @splitDuoData, fileName, AlgNames, opts);
% Ensure size of results does not include range over parameters
assert(all(size(results.obj_all) == [length(IndexConst.ObjectiveValueNames),...
                    opts.data_params.num_splits,numAlgs]));  
fprintf(1, '\n\n++++++++++++++++++++++++++++++++ Successful test for OVER_DEFAULT_PARAMS...\n\n');

% Ensure that can run parameters over best parameters;
fprintf(1, '\n\n++++++++++++++++++++++++++++++++ Running test for OVER_BEST_PARAMS...\n\n');
opts.cv_params.cross_validate = AlgConst.OVER_BEST_PARAMS;
[results, Xinfo, Yinfo, ModelInfo, params, bestParamIndices] = ...
    runBakeoff(@loadDuoData, @splitDuoData, fileName, AlgNames, opts);
% Ensure size of results includes all parameters
numParams = max(cellfun('length', params));
assert(all(size(results.obj) == ...
           [length(IndexConst.ObjectiveValueNames),2,numAlgs]));
assert(length(Xinfo) == opts.data_params.num_splits);
assert(length(Xinfo{1}) == IndexConst.NUM_VARS);
assert(length(Xinfo{1}{IndexConst.RECON_IND}) == numAlgs);
for a = 1:numAlgs
  assert(length(params{a}) >= bestParamIndices(a));
  assert(length(Xinfo{1}{IndexConst.RECON_IND}{a}) == length(params{a})); 
end
fprintf(1, '\n\n++++++++++++++++++++++++++++++++ Successful test for OVER_BEST_PARAMS...\n\n');


