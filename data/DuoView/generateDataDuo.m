function [X Y X1 X2 Y1 Y2] = generateDataDuo (n, c, t, opts)
% Now generates data descirbed in paper as well as shifted data
% Mag plays the role of sigma = opts.mag/100 for Gaussian noise
% And opts.ratio sparsifies B for B*H = X
% NOTE: Assumes all paths and global variables initialized

DEFAULTS = getDefaultDataParams();

if nargin < 3
    opts = DEFAULTS;
 else  
    opts = getOptions(opts, DEFAULTS);
end

B = randn(n, opts.rk);
nB = sqrt(sum(B.*B));
B = B ./ repmat(nB, n, 1);

W = randn(c, opts.rk);
nW = sqrt(sum(W.*W));
W = W ./ repmat(nW, c, 1);

Phi = rand(opts.rk, t);

X1 = B * Phi;
Y1 = W * Phi;

% MAKE SURE NOISE IS ADDED TO EACH ROW SO THAT EACH ROW HAS THAT OPTS.RATIO OF NOISE
% Add Gaussian noise to X is shifted and sparse noise to Y
if (opts.shifted == 1)
    sigma = opts.mag/100;
    X2 = sigma*randn(n,t);
  
	idx = randperm(c*t);
	Y2 = zeros(c, t);
	numSpY = floor(opts.ratio*c*t/100);
	Y2(idx(1:numSpY)) = rand(numSpY,1)*opts.mag;
       
else
	idx = randperm(n*t);
	X2 = zeros(n, t);
	numSpX = floor(opts.ratio*n*t/100);
	X2(idx(1:numSpX)) = rand(numSpX,1)*opts.mag;
	
	idx = randperm(c*t);
	Y2 = zeros(c, t);
	numSpY = floor(opts.ratio*c*t/100);
	Y2(idx(1:numSpY)) = rand(numSpY,1)*opts.mag;
end

fname = getDataNameDuo(n, c, t, opts);
X = X1 + X2;
Y = Y1 + Y2;
save(fname, 'X', 'Y', 'X1', 'X2', 'Y1', 'Y2');
