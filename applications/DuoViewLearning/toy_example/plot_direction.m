function [x1, x2, x3] = plot_direction(v, phi, opt_style)
% Makes the direction a unit direction

  if nargin < 3 || isempty(opt_style)
    persistent ind;
    styles = {'-.r*', '--cs', '-gd', '-.ko', '--m+',  '--y+'};
    if isempty(ind), ind = 1; end
    if ind > length(styles)
      ind = 1;
    end
    opt_style = styles{ind};
    ind = ind + 1;
  end
  
  v
  x1 = getPlotValue(v(1), phi);
  x2 = getPlotValue(v(2), phi);
  x3 = getPlotValue(v(3), phi);
  [x1 ; x2; x3]
  %  x1 = x1/norm(x1);  x2 = x2/norm(x2);  x3 = x3/norm(x3); 
  plot3(x1, x2, x3, opt_style);

  function x = getPlotValue(vi, phi)
    if abs(vi) < 1e-3
      x = zeros(size(phi));
    else
      x = phi/vi;
      x = x/norm(x)*8;
    end
  end

end
