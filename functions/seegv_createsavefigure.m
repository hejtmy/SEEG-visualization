function [ ] = seegv_createsavefigure( f, clims, clrmap, printResolution, outDir, figName)
%SEEGV_SLICE Summary of this function goes here
%   Detailed explanation goes here
% --- colorbar
    clrbar_axes = axes('visible','on', 'position', [0.97 0.2 0.01 0.6]);  % position
    clrbar_vals = [clims(1):1/1000 * diff(clims) : clims(2)]';
    clrbar_xAx = 2:10;
    clrbar_yAx = 1:size(clrbar_vals,1);
    clrbar_inds = cVals2cInds(repmat(clrbar_vals, [1, 10]), [clims(1), clims(2)], size(clrmap.brain, 1) + [1,size(clrmap.chnls, 1)]);
    clrbar_hndl = image(clrbar_xAx, clrbar_yAx, ones(size(clrbar_inds)));
    set(clrbar_hndl, 'CData', clrbar_inds, 'CDataMapping','direct');
    minInd = min(clrbar_yAx); maxInd = max(clrbar_yAx); 
    tickVals = round([minInd: (maxInd - minInd)/4 : maxInd]);                                                              % tick vals/labels
    tickName = cell(1,length(tickVals));
    for tick = 1:length(tickVals)
        tickName{tick} = num2str(clrbar_vals(tickVals(tick)), '%01.1f');
    end
    % set(clrbar_axes, 'XTick',[], 'YTick',tickVals, 'YTickLabel',tickName);
    set(gca, 'YDir', 'normal');
    set(gca, 'XTick', [], 'YTick', tickVals, 'YTickLabel', tickName);
    set(clrbar_axes, 'YAxisLocation','right');

    % --- output directory
    if ~exist('outDir', 'var')
        outDir = [pwd filesep 'slices'];
    end
    if ~exist(outDir, 'dir')
        mkdir(outDir);
    end  

    % --- figure name
    if ~exist('figName', 'var')
        figName = 'notNamed';
    end

    % --- save
    set(f, 'PaperPositionMode', 'auto');
    saveas(f, [outDir filesep figName '.fig']);
    
    if printResolution == 0
        print(f, '-dpng','-r0', [outDir filesep figName '.png']);
    else
        print(f, '-dpng','-r600', [outDir filesep figName '.png']);
    end
    close(f);
    display(['Figure: ' figName ' stored in: ' outDir]);
end