function [X,Y] = genGaussian(t,n,k,sigma)
% GENGAUSSIAN generates data X and Y
% if k = 1, then generate in -1, 1
% Since the Goldberg procedure is restrcted to k = 1,
% use genRevGaussian if k > 1
% All data assumptions are that X is stored as nxt, rather than txn!
% It'll be transposed later

% Goldberg
if k==1
    % Make sure that number of classes somewhat well balanced
    numTries = 6;
    maxDiff = 0.25;
    
    X = []; Y = [];
    i = 0;
    for i = 1:numTries
		% rank-r matrix
        r = n/2;
		X0 = randn(n,r)*randn(r, t);
        X0 = normalizeMatrix(X0');
        X0 = X0';
		X = X0 + sigma*randn(n,t);
		
		b = 0;
        W = 10*randn(1,n);
		y0 = W*X0 + b;
		r = rand(1,t);
		% P(1 | y0_ij)
		P = 1./(1+exp(-y0));
		Y = -1*ones(1,t);
		Y(r < P) = 1;
        
        X = X'; Y = Y';
        if abs(sum(Y))/t < maxDiff
            break;
        end
    end
    if i == numTries
        warning('genGaussian -> Could not evenly generate -1,1 classification, %d 1 and %d -1\n',sum(Y==1),sum(Y==-1));
    end    
else
    [X,Y,U] = genRevGaussian(t,n,k,sigma);
end


function [X,Y,U] = genRevGaussian(t,n,k,sigma)
% GENREVGAUSSIAN generates data X,Y and U
% if k = 1, then generate in -1, 1

% Generate random matrix of Y values Y1 = 1
Y = zeros(t,k);
if (k == 1)
	for i = 1:t
	    if (rand <= 0.5)
	        Y(i) = -1;
	    else
	        Y(i) = 1;
	    end    
	end	
else
    for i = 1:t
    	r = floor(k*rand)+1;
    	Y(i,r) = 1;
    end	
end

U = randn(k,n);
X = Y*U + sigma.*randn(t,n);

end


end
