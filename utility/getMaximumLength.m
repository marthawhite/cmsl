function [maxLen, maxIndex] = getMaximumLength(array)
% GETMAXIMUMLENGTH finds the maximum length of a subarray in array
% e.g.
%	>> A = {{{1},{2}}, {{3},{4},{5},{6}}, {{1}}}
%
% A = 
%
%   {1x2 cell}    {1x4 cell}    {1x1 cell}
%
% returns 4 with index 2
% Author: Martha White, University of Alberta, 2011
%
  maxLen = 0;
  maxIndex = 1;
  
  for i = 1:length(array)
    if length(array{i})> maxLen
      maxLen = length(array{i});
      maxIndex = i;
    end  
  end

end
