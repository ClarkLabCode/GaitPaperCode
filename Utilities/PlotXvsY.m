function PlotXvsY(x,y,varargin)
    %set up default values. all default values can be changed by varargin
    %by putting them in the command line like so
    %plotXvsY(...,'color','[0,0,1]');
    if(size(x,1) > 1)
        graphType = 'line';
    else
        graphType = 'scatter';
    end
    color = lines(size(y,2));
    error = [];
    lineStyle = '-';
    hStat = ishold;
    
    for ii = 1:2:length(varargin)
        eval([varargin{ii} '= varargin{' num2str(ii+1) '};']);
    end
    
    switch graphType
        case 'line'
            if isempty(error)
%                 plot(x,y,'lineStyle',lineStyle,'LineWidth',2);
                plot(x,y,'lineStyle',lineStyle,'LineWidth',1);
            else
                PlotErrorPatch(x,y,error,color);
            end
        case 'scatter'
            if isempty(error)
                scatter(x,y,50,color);
            else
                if ~hStat, hold on; end
                for c = 1:size(x,2)
                    scatter(x(:,c),y(:,c),50,color(c,:));
%                     errorbar(x(:,c),y(:,c),error(:,c),'color',color(c,:),'LineStyle','none','LineWidth',2);
                    errorbar(x(:,c),y(:,c),error(:,c),'color',color(c,:),'LineStyle','none','LineWidth',1);
                end
                if ~hStat, hold off; end
            end
        case 'bar'
            if ~hStat, hold on; end
            
            bar(x,y,'FaceColor',[0 0 1]);
            
            if ~isempty(error)
                errorbar(x,y,error,'LineStyle','none','LineWidth',2,'color',[0 0 0]);
            end
            
            if ~hStat, hold off; end
    end
end