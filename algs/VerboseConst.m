classdef VerboseConst
  properties (Constant)
		% All possibly verbose values, depending on how much
    % information is desired during algorithm runs
    % Higher numbers include lower values, i.e. detailed algorithm
    % information includes main script and basic algorithm info printed
    
    % No information printed out
    NONE = 0;
    % Information in runBaekoff or main script printed out    
    MAIN_SCRIPT = 1;
    MAIN_SCRIPT_COLOUR = 'black';
    % More detailed information in main script printed  
    MAIN_SCRIPT_DETAILED = 2;
    MAIN_SCRIPT_DETAILED_COLOUR = 'black';
    % Basic information, like termination tolerance,
    % printed by algorithm
    BASIC_ALG = 3;
    BASIC_ALG_COLOUR = 'blue';    
    % Detailed information, like minimization progress, 
    % printed by algorithm
    DETAILED_ALG = 4;
    DETAILED_ALG_COLOUR = 'blue';
    % Detailed information, like minimization progress, 
    % printed inside the generic solvers (adal, lbfgs, etc.)
    DETAILED_SOLVER = 5;
    DETAILED_SOLVER_COLOUR = 'red';
    % Basic information printing progress in cross validation
    BASIC_CV = 1;
    BASIC_CV_COLOUR = 'black';    
    % Detailed information about errors on folds
    DETAILED_CV = 2;
    DETAILED_CV_COLOUR = 'blue';
    
    WARNING_COLOUR = 'red';
  end
end