%TEST
clear;clc;

initialize_testing_paths;
addpath '../applications/DuoViewLearning';

filename = '../data/DuoView/all_faces_pos1_pos2_nx50_ny50_noiseoption1.mat';
if ~exist(filename)
  [X, Y, X1, Y1] = makeFaceData(50, 5, 1);
else  
 load filename;
end

X_noisy = X(:, 1:20);
Y_noisy = Y(:, 1:20);
X_clean = X1(:, 1:20); 
Y_clean = Y1(:, 1:20);
clearvars X X1 Y Y1;

global sigma_smooth_L11;
sigma_smooth_L11 = 1e-2;

[X_noisy,Y_noisy] = add_pixel_noise(X_noisy, X_clean, Y_noisy, Y_clean,...
                                    struct('percentage_ones', [0 0]));

reg_wgts = [0.1 1 10];
opts = [];
opts.verbose = VerboseConst.BASIC_ALG;
  
for r = reg_wgts
  opts.reg_wgt = r;
  filename = ['output_test_snr_approx_' num2str(r) '.mat'];
  if ~exist(filename)
    [X_recons, Y_recons] = convex_single_view(X_noisy, Y_noisy,opts);
    save(filename, 'X_recons', 'Y_recons');
  else
    load(filename);
  end  
  Z_recons = [X_recons;Y_recons];
  Z_noisy = [X_noisy;Y_noisy];
  Z_clean = [X_clean;Y_clean];
  
  % Test if the approximate SNR keeps the same order as true SNR
  Z_snr_approx = snr_approx(Z_recons, Z_noisy);
  Z_snr_approx2 = snr_approx(Z_recons);
  Z_recons_snr = snr(Z_clean,Z_recons);
  Z_noisy_snr = snr(Z_clean,Z_noisy);
  Z_recons_error = norm(Z_clean-Z_recons);
  Z_noisy_error = norm(Z_clean-Z_noisy);
  
  fprintf(1, '+++++++++++++++For reg_wgt %g:\n', opts.reg_wgt);
  Z_snr_approx
  Z_snr_approx2
  Z_recons_snr
  Z_noisy_snr
  Z_recons_error
  Z_noisy_error
end

