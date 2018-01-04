function [] = test()


x0 = [3, 2];    % Inital value
bl = [-2; 2];   % lower bounds
bu = [inf; 3];  % upper bounds


param.maxIter = 10;     % max number of iterations
param.maxLsIter = 20;   % max number of line search steps in each iteration
param.maxBdl = 10;      % max number of bundles to keep
param.maxFnCall = 100;  % max number of calling the function
param.tolCon = 1e-5;      % tolerance of constraint satisfaction
param.tolFun = 1e-5;   % final objective function accuracy parameter
param.verbose = 0;      % print intermediate steps (doesn't seem to work)

% param.dummy = 'dummy';

function [f, g] = toy_func(x)
    cen = [-1; 1];
    f = 0.5 * norm(x - cen)^2;
    g = x - cen;
end


% At present, only box constraints are implemented in the mex.
[x, f, iter, numCall, flag] = pbm(@toy_func, x0, bl, bu, param)





end
