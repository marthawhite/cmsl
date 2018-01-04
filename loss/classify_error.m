function err = classify_error(Yu, Output)
% assumes Yu is l x t
% is Yu in {-1,1}, converted to {0,1}

t = size(Yu, 2);

% Yu must be a 0/1 matrix.  It can be sparse for
% multiclass/multilabel
Yu(Yu == -1) = 0;
if any(Yu ~= 1 & Yu ~= 0)
    error('Ground truth labels must be 0/1');
end

err = 0;
for i = 1 : t
    soft_disc = Output(:, i);
    [max_disc, ori_idx] = max(soft_disc);
    
    idx = find(soft_disc >= max_disc);
    if length(idx) > 1
        err = err + 1;
    elseif Yu(ori_idx, i) == 0
        err = err + 1;
    end
end

err = err / t;
