function err = regression_error(Yu, Output)

[a, b] = size(Yu);
err = sum(sum((Yu - Output).^2)) / (a * b);
