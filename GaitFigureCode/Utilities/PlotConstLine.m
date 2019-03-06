function PlotConstLine(value,constDim)
        if nargin < 2
            constDim = 1;
        end

        referenceLineColor = [0.25,0.25,0.25];
        
        limitsX = xlim';
        limitsY = ylim';
        
        switch constDim
            case 1 % y value doesn't change
                % if the value is outside of the range, don't plot just
                % move the limits
                if (limitsY(1)>=value || limitsY(2)<=value)
                    sortedLimits = sort([limitsY; value]);
                    ylim([sortedLimits(1) sortedLimits(end)]);
                else
                    plot(limitsX,[value; value],'--','Color',referenceLineColor);
                end
            case 2 % x value doesn't change
                if (limitsX(1)>=value || limitsX(2)<=value)
                    sortedLimits = sort([limitsX; value]);
                    xlim([sortedLimits(1) sortedLimits(end)]);
                else
                    plot([value; value],limitsY,'--','Color',referenceLineColor);
                end
        end
    
end