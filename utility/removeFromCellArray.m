function newArray = removeFromCellArray(array, string)

ind = find(ismember(array, string));

newArray = array;
if ind > 0
  newArray(ind) = [];
end

end

