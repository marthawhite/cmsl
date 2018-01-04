function randInd = rand_int(range, numind)

lb = range(1);
ub = range(2);
randInd = floor((ub-lb+1)*rand(numind,1)+lb);

end
