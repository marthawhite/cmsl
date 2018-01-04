classdef AlgConst
  properties (Constant)
		% All specified algorithm names
    % 1. CSSL can also be specified as another name for
    %    'Multi-View', the global solver
    % 2. Single-View and Single-View-Row are solving the same thing
    %    for the smooth L11 loss, but single-view-row uses an algorithm
    %    developed specifically for the L11 loss in related work.
    % 3. Single-View-Col uses the same code as Single-View-Row (inexact_alm_pca)
    %    for a robust loss, BUT, it transposes the vectors to learn column-wise.
    % 4. The Staged alternator uses the same code as the three-way alternator
    %    over all the variables, but does not consider the L2 loss in the alterntion
    %    instead only optimizing B and Phi first; then, after convergence, it
    %    optimizes L2 with B and Phi fixed to obtain W.
    AllAlgNames = {'Alternator',  'Multi-View', 'Transductive-Matrix-Completion',...
                   'Single-View', 'Staged-Alt'};
    
    % Defaults for algorithm parameters
    DEFAULT_GAMMA = 1;
    DEFAULT_L2_WGT = 1;
    DEFAULT_NUM_BASIS = 10;
    DEFAULT_REG_WGT = 1e-3;
    
    % Ranges for algorithm parameters
    RANGE_GAMMA = [0.1 1];    
    RANGE_L2_WGT = [1e-3 1e-1 1];
    RANGE_NUM_BASIS = [10 50 100];
    RANGE_REG_WGT = [1e-3 1e-2 0.1];
    
    % Indices for CV parameter selection
    OVER_DEFAULT_PARAMS = 0;   % Returns results for default parameters 
    CV_PARAMS = 1;             % Cross validates over sets of params given by getParameters
    OVER_BEST_PARAMS = 2;      % Runs full experiment on all parameters, returns results
                               % struct for best parameters (but ModelInfo, etc, all runs)
    OVER_ALL_PARAMS = 3;       % Runs full experiments on all parameters, returns results
                               % struct for all parameters. Note that OVER_BEST_PARAMS can be
                               % used to all compute full results struct, as all need Xinfo
                               % Yinfo and ModelInfo is returned.
  end
end