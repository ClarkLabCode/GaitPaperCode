function ConfAxis(varargin)
    tickX = [];
    tickY = [];
    tickLabelX = [];
    tickLabelY = [];
    fTitle = [];
    figLeg = cell(0,1);
    labelX = [];
    labelY = [];
    fontSize = 20;

    for ii = 1:2:length(varargin)
        eval([varargin{ii} '= varargin{' num2str(ii+1) '};']);
    end
    
    if ~isempty(tickX)
        set(gca,'XTick',tickX)
    end
    
    if ~isempty(tickY)
        set(gca,'YTick',tickY)
    end
    
    % if length(XTickLabel) > length(XTick) matlab will ignore the extra
    % entries in XTickLabel
    if ~isempty(tickLabelX)
        set(gca,'XTickLabel',tickLabelX)
    end
    
    if ~isempty(tickLabelY)
        set(gca,'YTickLabel',tickLabelY)
    end
    
    if ~isempty(fTitle)
        title(fTitle,'fontSize',24);
    end
    
    if ~isempty(labelX)
        xlabel(labelX);
    end
    
    if ~isempty(labelY)
        ylabel(labelY);
    end
    
    if ~isempty(figLeg)
        legend(figLeg);
    end
    
    ax = gca;
    ax.YLabel.FontSize = fontSize;
    ax.XLabel.FontSize = fontSize;
    set(gca,'FontSize',fontSize)

%     if ~isempty(ax.Legend)
%         ax.Legend.FontSize = 20;
%     end
    set(gca,'LineWidth',2);
    set(gca,'box','off');
    
    %set(gca,'LooseInset',get(gca,'TightInset'))
end