
% Illustrates the usefulness of treating the two views
% separately, where only noise is added to the second view
% and the clean first view is useful in reconstructing the
% clean second view;
% Interestingly, MSL does not sacrifice accuracy on the clean
% X-View to reconstruct the Y-View well: it performs comparably
% on the x-view and is noticeably better on the y-view.

initialize_matlab;

global results_dir;

AlgNames = { 'Multi-View', 'Single-View'};

% Initialize parameters for experiment
% Make X have clean data and Y have the Gaussian
% noise already present in the RobutPCA data
opts = getDuoViewParams([0 -1]);
opts.verbose = VerboseConst.DETAILED_ALG;
%warning off;

opts.cv_params.cross_validate = AlgConst.OVER_DEFAULT_PARAMS;

opts.data_params.test = 0;
opts.data_params.tl = 200;
opts.data_params.num_splits = 10;

opts.problem_params.timeout = 5000;
opts.problem_params.recover = false;
opts.problem_params.override_params = true;
opts.problem_params.sigma = 1e-2;
opts.problem_params.reg_wgt = [0.05 0.1];
opts.problem_params.L2_wgt = 3;
opts.problem_params.gamma = 1;
opts.problem_params.num_basis = 20;

% Clear output file    
opts.output_file = [results_dir '/output_incfeatures'];
fileId = fopen(opts.output_file, 'w+');
fclose(fileId);

% Specify dataset sizes and regularization weights to iterate over
num_features = 20:20:100; %100:100:1000;

numDatasets = length(num_features);
numAlgs = length(AlgNames);
results.obj = zeros(2,numAlgs,numDatasets); 
results.rts = zeros(2,numAlgs,numDatasets);
results.err_x = zeros(2,numAlgs,numDatasets);
results.err_y = zeros(2,numAlgs,numDatasets);

% Set the parameters for the desired dataset to use
gopts = [];
gopts.rk = 10; gopts.ratio = 0.1; gopts.mag = 100; 
t = 5000; % Total points in dataset
    
allresults = [];
for d = 1:length(num_features)
  % Get data file
  f = num_features(d);
  fileName = get_data_name(f,f,t,gopts);
  if exist(fileName) == 0
    generate_data(f,f,t,gopts);
  else
    fprintf(1, 'Loading existing file %s\n', fileName);
  end  

  fileId = fopen(opts.output_file, 'a+');
  fprintf(fileId, '\n\n++++++++++++++++  Running with number of features = %g\n\n', f);
  fclose(fileId);  
   
	result = runBakeoff(@loadDuoData, @splitDuoData, ...
                       fileName, AlgNames, opts);
  allresults = [allresults {result}];
  
  % Normalize by number of features
  results.obj(:,:,d) = result.obj(1, :, :);
  results.rts(:,:,d) = result.rts(1, :, :);
  results.err_x(:,:,d) = result.err(2, :, :);
  results.err_y(:,:,d) = result.err(3, :, :);
end
LineNames = AlgNames;

% Save the results
%results_name = [results_dir '/results_incfeatures_' date];
results_name = [results_dir '/results_incfeatures'];
save(results_name, 'results', 'allresults', 'AlgNames',...
     'num_features','opts', 'gopts');

% Plot the training runtime results for increasing samples
plot_opts = [];
plot_opts.xtics = num_features;
plot_opts.xlabel = 'Number of Features';
plot_opts.fontSize = 24;
plot_opts.legend_options = {'Location', 'West', 'FontSize', plot_opts.fontSize-4};

plot_opts.title = 'Training Runtimes';
plot_opts.ylabel = 'Runtime (seconds)';
pdf_plot([results_dir '/figure_incfeatures_truntime'], squeeze(results.rts), LineNames, plot_opts);

% Plot the objective values for increasing samples
% Do not display legend, since same as for first graph
% plot_opts.legend_options = -1;
plot_opts.title = 'Objective Values';
plot_opts.ylabel = 'Objective Value';
pdf_plot([results_dir '/figure_incfeatures_objective'], squeeze(results.obj), LineNames, plot_opts);

% Plot the clean X error values for increasing samples
plot_opts.title = 'Clean X-Errors';
plot_opts.ylabel = 'L2 Reconstruction Error';
pdf_plot([results_dir '/figure_incfeatures_xerror'], squeeze(results.err_x), LineNames, plot_opts);

% Plot the clean Y error values for increasing samples
plot_opts.title = 'Clean Y-Errors';
plot_opts.ylabel = 'L2 Reconstruction Error';
pdf_plot([results_dir '/figure_incfeatures_yerror'], squeeze(results.err_y), LineNames, plot_opts);


