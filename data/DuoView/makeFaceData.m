function [X, Y, X1, Y1]  = makeFaceData (nx,ratio, zero_data)
% zero_data = 0 means add gaussian random noise
% zero_data = 1 means zero values in pictures
% zero_data = 2 means randomly

% Assumes that the directory Cropped_mat from the ExtendedYaleB dataset
% has been added to this directory. Get the datasets at:
% http://vision.ucsd.edu/extyaleb/CroppedYaleBZip/CroppedYale.zip
% Note that they are in pgm form; you will need to read them into
% Matlab and convert them to .mat
  
global data_dir;

if nargin < 3
    nx = 50;
    ratio = 5;
    zero_data = 1;
end
% downsample to nx*ny
ny = nx;
mag = 50;

root = [data_dir '/DuoView/ExtendedYaleB/Cropped_mat'];
root = '/Users/martha/Code/projects/sparse_coding/datasets/ExtendedYaleB/Cropped_mat'
rand('state', 1);   
randn('state', 1);

pids = [11:13 15 17:39];
poses = 1:2;
illums = 1:63;

t = 1000;    % number of pairs
n = nx*ny;
c = n;
X = zeros(n, t);
Y = zeros(n, t);
img_crop = [];

num_person = length(pids);
num_illum = length(illums);
num_pose = 2;

illu_id_X = 1;  % illumination id for view X
illu_id_Y = 2;  % illumination id for view Y

% loads faces as load([root '/yaleB11/crop_11_1_0.mat']);

% Sample the faces
for i = 1 : t

  person_id = pids(ceil(rand()*num_person));
  pose_id = (rand() > 0.5) + 1;
    
  fname = sprintf('%s/yaleB%d/crop_%d_%d_%d.mat', ...
          root, person_id, person_id, pose_id, illu_id_X);
  load(fname);  
  
  X(:, i) = downsample(double(img_crop(:)), n);
  fname = sprintf('%s/yaleB%d/crop_%d_%d_%d.mat', ...
          root, person_id, person_id, pose_id, illu_id_Y);
  load(fname);
  Y(:, i) = downsample(double(img_crop), nx, ny);
end

% the range of X and Y is [0, 255], so scale 
X = X / 255;
Y = Y / 255;
Y1 = Y;
X1 = X;

% Generate ratio percent of indices to be zero
total_num = n*t;
if zero_data == 1
    ratio = ratio/100;
    num_zeroed = floor(total_num*ratio);
    idx = randperm(total_num);

    X(idx(1:num_zeroed)) = 0;
    Y(idx(1:num_zeroed)) = 0;
elseif zero_data == 0
    ratio = ratio/100;
    X = X+ratio*randn(size(X));
    Y = Y+ratio*randn(size(Y));
else
    idx = randperm(total_num);
    X2 = zeros(n, t);
    numSpX = floor(ratio*total_num/100);
    X2(idx(1:numSpX)) = rand(numSpX,1)*mag;

    idx = randperm(total_num);
    Y2 = zeros(n, t);
    numSpY = floor(ratio*total_num/100);
    Y2(idx(1:numSpY)) = rand(numSpY,1)*mag;  
    
    X = X + X2;
    Y = Y + Y2;
end

save([data_dir '/DuoView/all_faces_pos1_pos2_nx' int2str(nx) ...
      '_ny' int2str(ny) '_noiseoption' int2str(zero_data) '.mat'],...
      'X', 'Y', 'X1', 'Y1');

end
