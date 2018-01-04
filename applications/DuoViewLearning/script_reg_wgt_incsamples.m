
% Tests the optimization performance of multi-view and alternator
% with a strict time limit.

initialize_matlab;

global results_dir;

AlgNames = { 'Multi-View', 'Alternator'};
msl_index = 1;
lsl_index = 2;

% Initialize parameters for experiment
opts = getDuoViewParams;
opts.result_params.error_fcns(1) = []; % Remove SNR as a measure

opts.verbose = VerboseConst.DETAILED_ALG;

opts.cv_params.cross_validate = AlgConst.OVER_DEFAULT_PARAMS;

opts.data_params.test = 0;
opts.data_params.num_splits = 2;

opts.problem_params.override_params = true;
opts.problem_params.sigma = 1e-2;

opts.solver_params.timeout = 200;
opts.solver_params.recover = false;
%opts.solver_params.init_with_singleview = 1;
%opts.solver_params.inner_solver = 'adal';

% Clear output file    
opts.output_file = [results_dir '/output_reg_wgts_incsample'];
fileId = fopen(opts.output_file, 'w+');
fclose(fileId);

% Specify dataset sizes and regularization weights to iterate over
tl = 10:10:50;
reg_wgts = [1e-3 5e-3 1e-2 0.05 0.1];
RegNames = {'\alpha=1e-3', '\alpha=5e-3', '\alpha=1e-2', '\alpha=5e-2','\alpha=1e-1'};

numDatasets = length(tl);
numAlgs = length(AlgNames);
numRegWgt = length(reg_wgts);
results.obj_all = zeros(opts.data_params.num_splits,numRegWgt,numDatasets); 
results.rts_all = zeros(opts.data_params.num_splits,numRegWgt,numDatasets);

% Set the parameters for the desired dataset to use
gopts = [];
gopts.rk = 10; gopts.ratio = 5; gopts.mag = 20; 
f = 100; % Number of features
t = 5000; % Total points in dataset
opts.problem_params.num_basis = f;  % Rank can be at most number of features

% Get data file
fileName = getDataNameDuo(f,f,t,gopts);
if exist(fileName) == 0
  generateDataDuo(f,f,t,gopts);
else
  fprintf(1, 'Loading existing file %s\n', fileName);
end
    
allresults = [];
for a = 1:length(reg_wgts)
  opts.problem_params.reg_wgt = reg_wgts(a);
  fileId = fopen(opts.output_file, 'a+');
  fprintf(fileId, '\n\n++++++++++++++++  Running increasing samples for regularizer weight %g\n\n', reg_wgts(a));
  fclose(fileId);  
  for i = 1:numDatasets
		opts.data_params.tl = tl(i);
   
		result = runBakeoff(@loadDuoData, @splitDuoData, ...
                         fileName, AlgNames, opts);
    allresults = [allresults {result}];
    
    % Save relative difference between LSL and MSL
    results.obj_all(:,a,i) = result.obj_all(1,:,lsl_index)./result.obj_all(1,:,msl_index);
    results.rts_all(:,a,i) = result.rts_all(1,:,lsl_index)./result.rts_all(1,:,msl_index);
  end 
end
results.obj = [mean(results.obj_all,1); std(results.obj_all,1)];
results.rts = [mean(results.rts_all,1); std(results.rts_all,1)];

% Plot the training runtime results for increasing samples
plot_opts = [];
plot_opts.xtics = tl;
plot_opts.xlabel = 'Number of Samples';
plot_opts.fontSize = 24;
plot_opts.legend_options = {'Location', 'West', 'FontSize', plot_opts.fontSize-4};

plot_opts.title = 'Training Runtime of LSL:MSL For Varying RegWgt';
plot_opts.ylabel = 'Runtime (seconds)';
pdf_plot([results_dir '/figure_lsl_vs_msl_training'], squeeze(results.rts), RegNames, plot_opts);

% Plot the objective values for increasing samples
% Do not display legend, since same as for first graph
%plot_opts.legend_options = -1;
plot_opts.title = 'Objective Values  of LSL:MSL For Varying RegWgt';
plot_opts.ylabel = 'Objective Value';
pdf_plot([results_dir '/figure_lsl_vs_msl_objective'], squeeze(results.obj), RegNames, plot_opts);




