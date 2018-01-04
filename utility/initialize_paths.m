function [] = initialize_paths ()

if exist('global_var', 'var')
    return;
end
global_var = 1;

% Directory variables
global data_dir;
global proj_dir;

% Function for printing progress
global col_print;

% Variable to determine if trying to use MEX; if USE_MEX = true
% and mex files are not used if possible, then warnings printed
global USE_MEX;
USE_MEX = false;

cd ..;
proj_dir = pwd;
cd utility;

% Add matlab solvers and loss
addpath(proj_dir);
addpath([proj_dir '/algs']);
addpath([proj_dir '/algs/recovery']);
addpath([proj_dir '/algs/competitors']);

addpath([proj_dir '/applications']);

addpath([proj_dir '/loss']);

addpath([proj_dir '/data/']);
addpath([proj_dir '/data/DuoView']);

addpath([proj_dir '/testing']);

addpath([proj_dir '/solver']);
addpath([proj_dir '/solver/adal']);
addpath([proj_dir '/solver/pbm']);
addpath([proj_dir '/solver/lbfgsb']);
addpath([proj_dir '/solver/uofa_optimizers']);
addpath([proj_dir '/solver/uofa_optimizers/nesterov']);
addpath([proj_dir '/solver/uofa_optimizers/svm']);
% For the solver variables, if lbfgs is not compiled
% then use lbfgs in uofa_optimizers
if exist(['lbfgsb.' mexext]) == 0
  addpath([proj_dir '/solver/uofa_optimizers/bfgs']);
  warning('Using slower, uncompiled LBFGS implementation.');
end
addpath([proj_dir '/solver/PROPACK']);

addpath([proj_dir '/utility']);
      
osname = getenv('OS');
if ~isempty(strfind(osname, 'Win'))
  col_print = @(varargin)(cprintf('Red', varargin{:}));
else
  col_print = @(varargin)(fprintf(1, varargin{:}));
end

% Specify the path to the data
data_dir = [proj_dir '/data'];


