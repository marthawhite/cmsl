initialize_matlab;

global data_dir;
global results_dir;

AlgNames = { 'Multi-View', 'Single-View', 'Alternator'};
ssl_index = 2;
msl_index = 1;
lsl_index = 3;

opts = [];
percentage_ones = [0 0.05];
opts = getDuoViewParams(percentage_ones);
opts.verbose = VerboseConst.DETAILED_SOLVER;

opts.cv_params.cross_validate = AlgConst.OVER_BEST_PARAMS;
opts.cv_params.comparison_index = 1; % y error against clean data

opts.data_params.test = 0;
opts.data_params.num_splits = 1;
opts.data_params.tl = 50;

opts.problem_params.override_params = true;
opts.problem_params.num_basis = 300;
opts.problem_params.sigma = 1e-3;
opts.problem_params.reg_wgt = [1 2 5 10];
opts.problem_params.gamma = 1;
opts.problem_params.L2_wgt = [2 5 10];

opts.solver_params.timeout = 1000;
opts.solver_params.recover = false;
opts.solver_params.inner_solver_lbfgs_params = [];
opts.solver_params.inner_solver_lbfgs_params.maxiter = 50;
opts.solver_params.solver_lbfgs_params = [];
opts.solver_params.solver_lbfgs_params.maxiter = 50;
opts.solver_params.maxiter_inner = 15;
%opts.solver.outer_solver_lbfgs_params = [];
%opts.solver.outer_solver_lbfgs_params.funTol = 1e-1;
%opts.solver.inner_solver_lbfgs_params =[];
%opts.solver.inner_solver_lbfgs_params.funTol = 1e-3;
%opts.solver.funTol = 1e-2;
%opts.solver.solver_lbfgs_params =[];
%opts.solver.solver_lbfgs_params.funTol = 1e-1;

nx = 100;
ny = nx;

fileName = [data_dir '/DuoView/all_faces_pos1_pos2_nx100_ny100_noiseoption1.mat'];
[results Xinfo Yinfo  ModelInfo params bestParamIndices] = ...
    runBakeoff(@loadDuoData, @splitDuoData, fileName, AlgNames, opts);
splitNum = 1;
index = 1;

% Plot True and Noisy Y image and then SSL and MSL below that
YImages = {Yinfo{splitNum}{IndexConst.CLEAN_IND}, Yinfo{splitNum}{IndexConst.NOISE_IND},...
           Yinfo{splitNum}{IndexConst.RECON_IND}{ssl_index}{bestParamIndices(ssl_index)},...
           Yinfo{splitNum}{IndexConst.RECON_IND}{lsl_index}{bestParamIndices(lsl_index)},...
           Yinfo{splitNum}{IndexConst.RECON_IND}{msl_index}{bestParamIndices(msl_index)}};
YLabels = {'Clean', 'Noisy','SSL', 'LSL', 'MSL'};
plotImage(YImages, YLabels, index, nx, ny, [results_dir '/figure_faces_' int2str(nx) '_' date]);

% Save information
saveLocation = [results_dir '/results_faces_' int2str(nx) '_' date];
save(saveLocation, 'results', 'AlgNames', 'opts', 'Xinfo', 'Yinfo', 'ModelInfo', 'params', 'bestParamIndices');

