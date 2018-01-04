function loss = l1_err(Z, Z_recons) 
% To test (X, X_recons, Y, Y_recons), pass Z-recons = [X_recons; Y_recons]
    loss = sum(sum(abs(Z - Z_recons)))/numel(Z);
end
