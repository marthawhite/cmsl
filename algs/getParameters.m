function AllParams = getParameters(algNames, getTypical)
% Returns all the parameter values to test over for the given algorithms
% Up to a max of iterating over 3 of the parameters. If algorithm has
% more than three paramaters to set, then it fixes one of the parameters.
% For example, for many sparse coding algorithms, the basis is fixed to the
% default to allow testing over L2_wgt, reg_wgt and gamma.
%
% If getTypical is set to 1, then only returns a "good" parameter option instead
%       of multiple options
% 
% For possible algorithms, see AlgConst.m
%
% Author: Martha White, University of Alberta, 2012

Algs = [];
numAlgs = length(algNames);
AllParams = [];
 
gammas = AlgConst.RANGE_GAMMA;  
L2_wgts = AlgConst.RANGE_L2_WGT;
numBasis = AlgConst.RANGE_NUM_BASIS;
reg_wgts = AlgConst.RANGE_REG_WGT;

if nargin == 2 && getTypical == 1
    L2_wgts = AlgConst.DEFAULT_L2_WGT;
    numBasis = AlgConst.DEFAULT_NUM_BASIS; 
    reg_wgts = AlgConst.DEFAULT_REG_WGT; 
    gammas = AlgConst.DEFAULT_GAMMA;
end

for ii = 1:numAlgs
    opts = [];  

    if strcmp(algNames{ii}, 'Alternator')       
        params = [];
        opts.num_basis = AlgConst.DEFAULT_NUM_BASIS;
        for a = reg_wgts
            for b = L2_wgts
                for g = gammas
                    opts.reg_wgt = a; opts.L2_wgt = b; opts.gamma = g;
                    params = [params {opts}];
                end
            end    
        end 
        AllParams = [AllParams {params}];   
     
    elseif strcmp(algNames{ii}, 'Transductive-Matrix-Completion')  
        params = [];
        for b = L2_wgts
            opts.L2_wgt = b;
            params = [params {opts}];
        end 
        AllParams = [AllParams {params}]; 

    elseif strcmp(algNames{ii}, 'Multi-View') || strcmp(algNames{ii}, 'CSSL')      
        params = [];
        opts.num_basis = AlgConst.DEFAULT_NUM_BASIS;
        for a = reg_wgts
            for L2_wgt = L2_wgts
                for g = gammas
                    opts.reg_wgt = a; opts.L2_wgt = L2_wgt; opts.gamma = g;
                    params = [params {opts}];
                end
            end   
        end 
        AllParams = [AllParams {params}]; 

    elseif strcmp(algNames{ii}, 'Single-View')
        params = [];
        for a = reg_wgts
          for g = gammas
            opts = [];
            opts.reg_wgt = a; opts.gamma = g;
            params = [params {opts}];
          end
        end 
        AllParams = [AllParams {params}];   
     
	 elseif strcmp(algNames{ii}, 'Staged-Alt')       
        params = [];
        for a = reg_wgts
            for n = numBasis
                for g = gammas
                    opts.L2_wgt = 0; opts.reg_wgt = a; opts.num_basis = n; opts.gamma = g;
                    params = [params {opts}];
                end
            end
        end 
        AllParams = [AllParams {params}]; 
        
    else
        error(['getAlgs -> Cannot currently handle ' algNames{ii} ' must use ' ...
               'predefined algorithms']);
        AllAlgNames
    end

end

