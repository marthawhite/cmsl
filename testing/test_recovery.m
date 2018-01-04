%TEST
% Driver script to test recovery algorithms
warning on;

initialize_testing_paths
addpath ../algs/recovery/old_recover;

rand('state', 1);
randn('state', 1);

% the matrix Zhat is (n+c)*t
t = 100;
n = 100;
c = 100;
r = 20;   % rank of Zhat (for method 3 of generating Zhat)

gamma = 1;

% Generate the matrix Zhat
m = n+c;

gen_data_method = 1;

for gen_data_method = 1:3
  fprintf(1, '***********************************************************\n');
  fprintf(1, 'Generating data with method %u\n', gen_data_method);
  switch gen_data_method
  case 1
	  Z = rand(m, t);
	case 2
	  Z = randn(m,r) * randn(r,t);
	case 3
	  Phi = randn(r, t);
	  B = randn(n, r);
	  W = randn(c, r);
	  for i = 1 : r
	    B(:,i) = B(:,i) / norm(B(:,i));
	    W(:,i) = W(:,i) / (norm(W(:,i)) / gamma);
    end
    U = [B; W];
	  Z = U*Phi;
  end

  % Invoke the recovery routine
  diff_thrd = 1e-2;     % exiting threshold of relative error in Frobenius norm
  dims = [n; c; t];

  fprintf(1, '---------------- Using cone recovery on data %u\n', gen_data_method);
  start_time = cputime;
  [Uhat1 Phihat1] = cone_recover(Z, struct('gamma', gamma, 'verbose', 0), dims, diff_thrd);
  fprintf(1, 'Cone recover took %g seconds\n', cputime - start_time);
  fprintf(1, 'Cone recover error on Z decomposition: %g\n', norm(Z - Uhat1*Phihat1));
  if gen_data_method == 3
    fprintf(1, 'Cone recover rank: %u, true rank: %u\n', rank(Phihat1), r);
    %fprintf(1, 'Cone recover error on representation: %g\n', norm(Phi - Phihat1));
  end
  
  fprintf('\n\n---------------- Using old recovery on data %u\n', gen_data_method);
  start_time = cputime;
  [Uhat2 Phihat2] = old_recover(Z,dims,gamma);
  fprintf(1, 'Old recover took %g seconds\n', cputime - start_time);
  fprintf(1, 'Old recover error on Z decomposition: %g\n', norm(Z - Uhat2*Phihat2));
  if gen_data_method == 3
    fprintf(1, 'Old recover rank: %u, true rank: %u\n', rank(Phihat2), r);
    %fprintf(1, 'Old recover error on basis: %g\n', norm(U - Uhat2));
    %fprintf(1, 'Old recover error on representation: %g\n', norm(Phi - Phihat2));
  end
  
end

