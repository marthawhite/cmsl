% Constant return errors for both uofa lbfgs implementation
% and mex implementation
classdef BFGS_ERRORS
  properties (Constant)

    % Possible error values for flag
    % Made large because trying not to overlap
    % with lbfgsb which has its own flags
    SUCCESS = 1;    % No errors
    MAXITER = 5;    % Reached max iterations in main loop
    TIMEOUT = 101;    % Reached max time in main loop
    CONSTRAINT = 102; % Error in enforcing constraints
    BACKTRACK_ITER = 103;  % Max iters in backtrack
    BACKTRACK_ERROR = 104;  % Error in backtrack

    % Other MEX implementation error names
    SMALL_REL_CHANGE = 2;
    SMALL_PROJ_GRADIENT = 3;
    MAX_FCN_CALL = 4;
    INPUT_ERROR = -1;
    ABNORMAL_TERMINATE = -2;
  end
end