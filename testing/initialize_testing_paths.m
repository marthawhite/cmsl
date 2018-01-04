function [] = initialize_testing_paths ()

	cd ../utility
	initialize_paths
	cd ../testing
  
  global sigma_smooth_L11;
  sigma_smooth_L11 = 1e-2;

end

