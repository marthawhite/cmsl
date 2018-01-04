function plotImage(X, Xlabel, index, nx, ny, pdf_name, width, blankSpots)
  
  % If X is a cell array, then it contains the X's to be plotted as subplot
  if (iscell(X))
    if nargin < 7, width = 3; end    
    if nargin < 8, blankSpots = []; end    
    numPlots = ceil((length(X)+length(blankSpots))/width);
    i = 1;
    ind = 1;
    figure;
    for i = 1:length(X)+length(blankSpots)
      if ~isempty(blankSpots) && ~isempty(find(blankSpots == i)), continue; end
      subplot(numPlots,width,i);
      plotPic(X{ind}(:,index));
      title(Xlabel{ind},'FontWeight','bold', 'FontSize',16);
      ind = ind+1;
    end
  else
    % Plot single image
    figure; 
    plotPic(X(:, index));
  end
  
  if nargin >= 6
    hold on;
    name = [ pdf_name '.eps'];
    saveas(gcf, name, 'psc2');
    hold off;
    fprintf(1, 'Outputted figure to %s\n', name);    
  end
  
  function plotPic(Image)
    imagesc(uint8(round(255*reshape(Image,[nx,ny]))), [0, 255]); colormap(gray);
  end    
end
