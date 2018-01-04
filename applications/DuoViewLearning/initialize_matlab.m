function [] = initialize_matlab ()


if exist('global_var', 'var')
    return;
end
global_var = 1;

cd ../../utility
initialize_paths;
cd ../applications/DuoViewLearning/
  
global proj_dir;
global results_dir;
results_dir = [proj_dir '/applications/DuoViewLearning/results'];

