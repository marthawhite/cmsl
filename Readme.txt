*********************************************************************************
Folder Overview:

algs/:  the implementations of our algorithm and competitors. This folder also includes the recovery algorithm for the global multi-view and single-view solvers. Note that the global multi-view algorithm can also use the new boosting approach, which does not require a separate recovery procedure. That algorithm will be included in the next release.

applications/: specific scripts and settings to run experiments for each of the settings (multi-view and semi-supervised learning). The main file for running all scripts is runBakeoff.m, which runs all specified algorithms against each other for given settings. Currently only contains scripts for multi-view; scripts for semi-supervised learning are coming soon.

data/: the datasets for each of the different settings (e.g. semi-supervised learning and multi-view learning), which can include synthetic dataset generators.

loss/: Some loss functions for classification, input reconstruction and regularizers (like the trace norm).

solver/: the main solvers, including nesterov, adal, lbfgs, pbm and boosting solvers. The folder uofa_optimizers contains optimizers we have written. The folder lbfgsb contains a mex-file implementation of lbfgs by another group. The default is to use the UofA implementation.

testing/: some testing routines. Run test_all.sh to run all the tests.

utility/: some utility routines, such as for plotting.

*********************************************************************************
Getting Started:

You can easily get started by running one of the scripts in applications/DuoViewLearning, which should work out-of-the-box. For example, run script_incfeatures to compare the convex multi-view, the alternating multi-view and the convex single-view algorithms, for an increasing number of features. To change the number of samples, parameters, etc. used, go into script_incfeatures to see where they are set. To change the number of training samples, for example, change opts.tl to whatever amount you might want.

To get better performance, you can compile some mex files in solver/lbfgsb and solver/pbm.
I find that on a Mac, typing 'make all' in solver/lbfgsb works and on a Linux, running 'make' in Matlab seems to work. If you want to downsample the face data in data/DuoView, then you have to type 'mex downsample.cpp' in Matlab.
*********************************************************************************
More Detailed Information:

+++ Passing in many parameter options:
There is some general commonality amongst the code. Generally, each file takes in some parameters from an 'opts' variables. This opts variable has DEFAULT values in each file, and whatever is specified in the passed opts overrides the default value. For example, in the solver 'solvers/uofaoptimizers/nesterov/solve_Nesterov_generic.m':

  DEFAULTS.init_L = 1 / size(X0,2);  % estimate of Lipschitz constant, adjusted during opt
  DEFAULTS.reg_wgt = 1;  % weight on regularizer on Phi
  DEFAULTS.maxiter = 100; % previously 200, changed for time issues
  DEFAULTS.funTol = 1e-6;
  DEFAULTS.verbose = 0;

Passing in opts.reg_wgt = 0.1 will override that default, i.e.
solve_Nesterov_generic(loss_fun, reg_fun, prox_op, X0, struct('reg_wgt', 0.1))
causes the following opts to be used in solver_Nesterov_generic:
  opts.init_L = 1 / size(X0,2); 
  opts.reg_wgt = 0.1;
  opts.maxiter = 100;
  opts.funTol = 1e-6;
  opts.verbose = 0;

To determine which parameters can be passed to each algorithm, simply go to the file and look at the DEFAULTS that are specified.

+++ Explanation of algorithms
getAlgs and getParameters returns all the possible algorithms and their corresponding range or default parameters. For cross validation, a range of parameters is returned by getParameters; to change the range, modify getParameters.

Note that only basic parameters for the algorithms are set in getParameters (e.g. regularizer weight), the losses are set in runBakeoff and the remaining values take on the defaults in each algorithm file. To set maxiter or other 'opts' in each algorithm besides the ones returned by getParameters, use 'addToAllOptions'. This function takes the current array of parameter options (which will be cross validated over) and adds all the fields in  the second options struct passed (see addToAllOptions for details).

+++ Detailed Explanation of one Multi-View Script

To understand how to use some of the functionality, this section explains applications/DuoViewLearning/script_incfeatures.m. 
1. initialize_matlab calls initialize_paths.m in /utility to add all the paths to solvers, etc. It also sets the global variables, such as the location of the data.
2. AlgNames specifies which algorithms you want run. See AlgConst.m for a list of all possible algorithms.
3. getDuoViewParams initializes many parameters, including errors and losses, for the duo-view subspace learning setting. For example, its sets the add_noise function in opts to be 'add_pixel_noise': this function takes the percentage of ones and changes the noise. If the percentage_ones = [0 5], then it returns a clean X dataset for learning and a Y data with 5% pixel noise. If the percentage_ones = [5 -1], then X has pixel noise and Y has the original noise given in the dataset (in this case, Gaussian noise). Other settings are that the errors reported are the SNR, Clean-X-Error and Clean-Y-Error. 
4. The opts variable specifies the parameters for the experiment. The key parts are opts.cv_params, opts.data_params, opts.problem_params and opts.results_params. See runBakeoff.m for all possible experiment settings.
5. On line 37, we specify the different datasets we will run over: ones with an increasing number of features. We will record the objective value, runtimes and error between the reconstructed X and Y and their clean counterparts. The variable 'gopts' specifies the parameters for the data generation (e.g. gopts.ratio gives amount of noise added). 
6. script_incfeatures.m then runs through the features sizes given in num_features. 
7. Afterwards, the results are saved and plotted. Plotting is facilitated with pdf_plot, which prints the plot to a pdf if a Matlab desktop is not being used (i.e. matlab -nodisplay is run). See pdf_plot for options to set the fontSize, legend, etc.

+++ Detailed Explanation of Datasets

Synthetic data is generated on the fly. So, if you set certain data options (number of samples, noise, etc.) in you script (e.g. script_incfeatures.m), and a dataset with those options does not exist, then it will be generated. Once that dataset exists, it is not generated again. To regenerate new data, simply remove that dataset.

*********************************************************************************
Comments:

1. The local alternator was also implemented, with some optimizations, by Lee et al. (2009) for the paper 'Exponential family sparse coding with applications to self-taught learning'. We did not find a significant difference between our alternator and theirs (which used the staged approach, i.e. alternator_multi_view.m with opts.L2_wgt = 0). Please check out their code for your own comparison.

2. Future releases may contain an implementation of Bayesian canonical correlation analysis; preliminary experiments do not indicate good performance, so for now, it is left out of the competitors.

3. There have been single-view solvers designed for the L11 loss, by the Perception and Decision Laboratory, University of Illinois, Urbana-Champaign, Microsoft Research Asia, Beijing. There is an implementation given by Minming Chen, October 2009. In our experiments, these implementations did not outperform our single-view solver (and are specific to robust PCA), so we do not include them.
