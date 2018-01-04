function [] = pdf_plot(pdf_name,results,LineNames,opts)
% Outputs the plot to a pdf file
% results is a 2 X numAlgs X numPoints maxtrix 
% where the first row contains
% the means and second row the confidence interval
% opts.xlabel = 'Dataset Size' by default
% opts.ylabel = 'Accuracy' by default
% opts.datasetName = "Synthetic Gaussian" by default
% opts.display_legend = 1 by default

DEFAULTS.xlabel = 'Dataset Size';
DEFAULTS.ylabel = 'Misclassification Error';
DEFAULTS.xtics = [];
DEFAULTS.xvals = [];
DEFAULTS.title = [];
DEFAULTS.legend_options = [];
DEFAULTS.ignore_alg_ind = [];
DEFAULTS.fontSize = 28;
DEFAULTS.lineWidth = 0.8;
DEFAULTS.lineStyles = [];
DEFAULTS.annotate = 0;

if nargin < 4
  opts = DEFAULTS;
else  
  opts = getOptions(opts, DEFAULTS);
end

% Get subset of desired indices
desired_algs = [];
if ~isempty(opts.ignore_alg_ind)
  for i = 1:length(AlgNames)
    if isempty(find(opts.ignore_alg_ind == i)) 
      desired_algs = [desired_algs i];
    end    
  end
else
  desired_algs = 1:length(LineNames);
end

% Keep styles the same for each algorithm
% Only used for extra algorithms
if isempty(opts.lineStyles)
  [opts.lineStyles, LineNames] = getLineStyles(LineNames);
end

figure
if ~isempty(opts.xtics) && length(opts.xtics) == size(results,3)
  if ~iscell(opts.xtics)
    opts.xtics = vec2strcell(opts.xtics);
  end    
  set(gca,'XTick',1:length(opts.xtics));
  set(gca,'XTickLabel',opts.xtics);
end    
set(gca, 'FontSize', opts.fontSize);
xlabel(opts.xlabel,'FontSize',opts.fontSize);
yhandle = ylabel(opts.ylabel,'FontSize',opts.fontSize);
if ~isempty(opts.title)
  % Align title with y-axis label
  titlehandle = title(opts.title,'FontWeight','bold','FontSize',opts.fontSize);
  postitle = get(titlehandle,'Position');
  posy = get(yhandle,'Position');
  set(titlehandle,'Position',postitle);
end


hold
if ~isempty(opts.xvals)
  for i = desired_algs
    errorbar(opts.xvals, results(1,i,:), results(2,i,:), opts.lineStyles{i}, ...
             'LineWidth',opts.lineWidth, 'MarkerSize',10);
  end    
else    
  for i = desired_algs
    errorbar(results(1,i,:), results(2,i,:), opts.lineStyles{i}, ...
             'LineWidth',opts.lineWidth, 'MarkerSize',10);
  end
end

newLineNames = {LineNames{desired_algs}};
if opts.annotate == 1
  for i = desired_algs
    text(opts.xvals(end), results(1,i,end), ['\leftarrow' LineNames{i}],...
         'HorizontalAlignment','left')
  end
elseif isempty(opts.legend_options) 
  legend(newLineNames);
elseif iscell(opts.legend_options)
  ll = legend(newLineNames);
  for ii = 1:2:length(opts.legend_options)-1
    set(ll, opts.legend_options{ii}, opts.legend_options{ii+1});
  end  
end

name = [pdf_name '.eps'];
saveas(gcf, name, 'psc2');

function str_array = vec2strcell(vector)
  str_array = [];
  for i = 1:length(vector)
    str_array = [str_array {num2str(vector(i))}];
  end
end

function [STYLES, LineNames] = getLineStyles(LineNames)
% If LineNames is algorithm names, then always assigns
% an algorithm to the same line colour
SYMBOLS = {'-.r*', '-gd', '-.ko', '--m^',  ':cs', 'y+'};
extraInd = 1;
STYLES = [];
numAlgs = length(LineNames);
for a = 1:numAlgs
  algName = LineNames{a};
	switch algName
	  case 'Staged-Alt'
	    STYLES = [STYLES {SYMBOLS{1}}];
	  case {'Alternator', 'LSL'}
	    STYLES = [STYLES {SYMBOLS{2}}];  
      LineNames{a} = 'LSL';
	  case {'Single-View', 'SSL'}
	    STYLES = [STYLES {SYMBOLS{3}}];
      LineNames{a} = 'SSL';      
    case {'Multi-View', 'MSL', 'CSSL'}
	    STYLES = [STYLES {SYMBOLS{4}}];
      LineNames{a} = 'MSL';
	  case {'Single-View-Row', 'SSL-Row'}
	    STYLES = [STYLES {SYMBOLS{5}}];      
	  case {'Single-View-Col', 'SSL-Col'}
	    STYLES = [STYLES {SYMBOLS{6}}]; 
	  case 'FPCB'
	    STYLES = [STYLES {SYMBOLS{7}}];  
	  otherwise
	    STYLES = [STYLES {SYMBOLS{extraInd}}];
      extraInd = extraInd + 1;
 end
end  
end

end
 
