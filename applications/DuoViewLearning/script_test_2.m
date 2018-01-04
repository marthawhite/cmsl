initialize_matlab;

global data_dir;
global results_dir;

AlgNames = { 'Multi-View', 'Single-View', 'Alternator'};
AlgNames = { 'Multi-View', 'Single-View'};
ssl_index = 2;
msl_index = 1;
lsl_index = 1;

opts = [];
percentage_ones = [0 0.05];
opts = getDuoViewParams(percentage_ones);
opts.verbose = VerboseConst.DETAILED_SOLVER;

opts.cv_params.cross_validate = AlgConst.OVER_BEST_PARAMS;
opts.cv_params.comparison_index = 3; % y error against clean data
opts.cv_params.max_samples = 5; 
opts.cv_params.num_folds = 5; 
opts.cv_params.verbose = VerboseConst.DETAILED_CV; 
opts.cv_params.timeout = 50; 

opts.data_params.test = 0;
opts.data_params.num_splits = 1;
opts.data_params.tl = 100;

opts.problem_params.override_params = true;
opts.problem_params.num_basis = 200;
opts.problem_params.sigma = 1e-3;
opts.problem_params.reg_wgt = [0.1 0.5 1];
opts.problem_params.gamma = 1;
opts.problem_params.L2_wgt = [1 2 5];

opts.solver_params.timeout = 500;
opts.solver_params.recover = false;
opts.solver_params.inner_solver_lbfgs_params = [];
opts.solver_params.inner_solver_lbfgs_params.maxiter = 50;
opts.solver_params.solver_lbfgs_params = [];
opts.solver_params.solver_lbfgs_params.maxiter = 50;
opts.solver_params.maxiter_inner = 15;

nx = 50;
ny = nx;
%opts.solver_params.latent_dim = 30;

%rand('state', 3);
%randn('state', 3);

fileName = '/Users/martha/Code/projects/sparse_coding/releases/CSSL-2.0/data/DuoView/all_faces_pos1_pos2_nx50_ny50_noiseoption1.mat';
[results, Xinfo, Yinfo, ~, ~, bestParamIndices] = ...
    runBakeoff(@loadDuoData, @splitDuoData, fileName, AlgNames, opts);
splitNum = 1;
index = 1;

% Plot True and Noisy Y image and then SSL and MSL below that
YImages = {Yinfo{splitNum}{IndexConst.CLEAN_IND}, Yinfo{splitNum}{IndexConst.NOISE_IND},...
           Yinfo{splitNum}{IndexConst.RECON_IND}{ssl_index}{bestParamIndices(ssl_index)},...
           Yinfo{splitNum}{IndexConst.RECON_IND}{lsl_index}{bestParamIndices(lsl_index)},...
           Yinfo{splitNum}{IndexConst.RECON_IND}{msl_index}{bestParamIndices(msl_index)}};
YLabels = {'Clean', 'Noisy','SSL', 'LSL', 'MSL'};
plotImage(YImages, YLabels, index, nx, ny, [results_dir '/figure_faces_test']);


