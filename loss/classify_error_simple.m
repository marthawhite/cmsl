function err = classify_error_simple(Yu, Output)
% assumes Yu is 1 x t in {-1,1}

t = size(Yu, 2);
assert(size(Yu, 1) == 1);
Yu = sign(Yu);

% Yu must be a -1/1 vector
if any(Yu ~= -1 & Yu ~= 1)
    error('Ground truth labels must be -1/1');
end

err = sum(min(abs(Yu-Output),1))/t;

end

