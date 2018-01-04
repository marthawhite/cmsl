%TEST
clear;clc;

initialize_testing_paths;

% remove TMC as it only works for semi-supervised learning
TOL = 1e-3;
AlgNames = AlgConst.AllAlgNames;
AlgNames = removeFromCellArray(AlgNames, 'Transductive-Matrix-Completion');
Algs =  getAlgs(AlgNames);
Parameters = getParameters(AlgNames, 1);
for alg = 1:length(Algs)
  Parameters{alg}{1}.verbose = 1;
  err = factorization_test(Algs{alg}, Parameters{alg}{1});
  if abs(err) > TOL
    cprintf('red', ['Factorization test failed for %s, obtained error = %g'...
            ', > TOL = %g\n'], AlgNames{alg}, err, TOL);
  end  
end

