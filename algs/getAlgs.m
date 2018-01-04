function Algs = getAlgs(algNames)
% Returns all the function pointers in the format
%
%   [X_learned, Y_learned, pobj, runtime, B, W, Phi] = fcn(X, Yl, opts) 
%
%   where X = [Xl Xu], Yl is labeled data (or another view) and opts
%   are algorithms options. Note that for transductive algorithms,
%   B, W and Phi are not returned.
%
% For possible algorithms, see AlgConst.m
%
% Author: Martha White, University of Alberta, 2011

Algs = [];
numAlgs = length(algNames);

for ii = 1:numAlgs
    
  if strcmp(algNames{ii}, 'Alternator')       
    fcn = @alternator_multi_view;
    Algs = [Algs {fcn}]; 
    
  elseif strcmp(algNames{ii}, 'Transductive-Matrix-Completion')       
    fcn = @transductive_matrix_completion;
    Algs = [Algs {fcn}]; 

  elseif strcmp(algNames{ii}, 'Multi-View') || strcmp(algNames{ii}, 'CSSL')   
    fcn = @convex_multi_view;
    Algs = [Algs {fcn}]; 

  elseif strcmp(algNames{ii}, 'Single-View')
    fcn = @convex_single_view;
    Algs = [Algs {fcn}];        
 
  % Staged alternator simply alternator with beta = 0
  % Where opts.L2 optimized after 'opts.L1 + reg_wgt reg_loss' optimized.
	elseif strcmp(algNames{ii}, 'Staged-Alt')       
    fcn = @alternator_multi_view;
    Algs = [Algs {fcn}]; 
    
  else
    error(['getAlgs -> Cannot currently handle ' algNames{ii} ' must use ' ...
           'predefined algorithms']);
    AllAlgNames
 end

end

