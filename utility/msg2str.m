function msg_str = msg2str(msg, solver_name, opts, added_string)
  
  if nargin < 4, added_string = ''; end
  
  if strcmp(solver_name, 'pbm') == 1
    switch msg
      case 0
        msg_str = 'The problem has been solved.';
      case 1
        msg_str = 'Number of calls of function = maxFnCall.';
      case 2
        msg_str = 'Number of iterations = maxIter.';
      case 3
        msg_str = 'Invalid input parameters.';
      case 4
        msg_str = 'Not enough working space.';
      case 5
        msg_str = 'Failure in quadratic program';
      case 6
        msg_str = 'The starting point is not feasible.';
      case 7
        msg_str = 'Failure in attaining the demanded accuracy.';
      case 8
        msg_str = 'Failure in function or subgradient calculations.';
      otherwise
        msg_str = ['Unknown exit code: ' num2str(msg) 'in pbm'];
    end
  end
    
  if strcmp(solver_name, 'lbfgs') == 1
    switch msg
      case BFGS_ERRORS.SUCCESS
        msg_str = 'The problem has been solved.';
      case BFGS_ERRORS.SMALL_REL_CHANGE
        msg_str = 'The relative change of obj value < relative change.';
      case BFGS_ERRORS.SMALL_PROJ_GRADIENT
        msg_str = 'The norm of projected gradient < tolPG.';
      case BFGS_ERRORS.MAX_FCN_CALL
        msg_str = 'Number of function calls > maxFnCall.';
      case BFGS_ERRORS.MAXITER
        if nargin >= 3 && isfield(opts, 'maxiter')
          msg_str = ['Number of iterations > maxiter = ' num2str(opts.maxiter) added_string];
        else  
          msg_str = 'Number of iterations > maxIter\n';
        end    
      case BFGS_ERRORS.INPUT_ERROR
        msg_str = 'The routine has detected an error in the input parameters.';
      case BFGS_ERRORS.ABNORMAL_TERMINATE
        msg_str = 'Terminated abnormally and was unable to satisfy the convergence criteria.';
      case BFGS_ERRORS.TIMEOUT
        if nargin >= 3 && isfield(opts, 'timeout')
          msg_str = ['Reached maximum alloted time for optimization = '...
                     num2str(opts.timeout) added_string];
        else  
          msg_str = 'Reached maximum alloted time for optimization\n';
        end       
      case BFGS_ERRORS.BACKTRACK_ITER
        if nargin >= 3 && isfield(opts, 'backtrack_maxiter')
          msg_str = ['Maximum number of iterations in backtrack search = '...
                     num2str(opts.backtrack_maxiter) '\n'];
        else
          msg_str = 'Maximum number of iterations in backtrack search\n';
        end
      case BFGS_ERRORS.BACKTRACK_ERROR
          msg_str = 'Error in backtrack search';
      case BFGS_ERRORS.CONSTRAINT
        msg_str = 'Error enforcing constraints';           
      otherwise
        msg_str = ['Unknown exit code: ' num2str(msg) ' in lbfgs'];
    end
  end
end