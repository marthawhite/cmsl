function [Xi_noisy, Yi_noisy, Xi_clean, Yi_clean, Xtest, Ytest, Xtest_clean, Ytest_clean] = ...
      splitDuoData(X_noisy, Y_noisy, X_clean, Y_clean, idx, opts)  

  num_ex = size(X_noisy,2);
  idx = randperm(num_ex);    
  idxtrain = idx(1:opts.data_params.tl);
  opts.data_params.test = min(opts.data_params.test, num_ex);
  idxtest = idx((end - opts.data_params.test+1):end);
  
  Xi_noisy = X_noisy(:,idxtrain);
  Yi_noisy = Y_noisy(:,idxtrain);
  Xi_clean = X_clean(:,idxtrain);
  Yi_clean = Y_clean(:,idxtrain);
  
  if opts.problem_params.recover
    Xtest = X_noisy(:, idxtest);
    Ytest = Y_noisy(:,idxtest); 
    Xtest_clean = X_clean(:, idxtest);
    Ytest_clean = Y_clean(:,idxtest);   
  else
    Xtest = Xi_noisy; Ytest = Yi_noisy; Xtest_clean = Xi_clean; Ytest_clean = Yi_clean;
  end  

end