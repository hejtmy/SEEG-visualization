function [ output_args ] = seegv_slicessagital(mni_coors, dim )
%SEEGV_SLICESSAGITAL Summary of this function goes here
%   Detailed explanation goes here

sliceStep = 1.5;            % in [mm]
n_maxSlices = 120;
dim = 1;

min_mni = min(mni_coors(:, dim));
max_mni = max(mni_coors(:, dim));
i_sliceStep = closestval(xi, sliceStep) - closestval(xi, 0.0);
i_slices = closestval(xi, min_mni):i_sliceStep:closestval(xi, max_mni);    
n_slices = length(i_slices);
while n_slices > n_maxSlices
    sliceStep = sliceStep + voxSize_new;        % increase slice step
    i_sliceStep = closestval(xi, sliceStep) - closestval(xi, 0.0);
    i_slices = closestval(xi, min_mni):i_sliceStep:closestval(xi, max_mni);    
    n_slices = length(i_slices);
end
% --- figure
f = figure('visible','on');
set(f, 'Position', plotInfo.figurePosition);
%set(f, 'Position', [1 -479 2880 1472]);
set(f, 'Colormap', clrmap.fig);

nRows = floor(sqrt(n_slices/2));    % good coverage when: nCols = 2 * nRows
nCols = ceil(2*sqrt(n_slices/2));
if n_slices > nRows * nCols, nRows = nRows+1; end
assert(n_slices <= nRows * nCols);
marg_h = [0.05 0.05];
marg_w = [0.04 0.04];
gap = [0.005, 0.005];

% --- text
if isfield(plotInfo, 'text')
    tx = axes('visible','off', 'position',[0 0 1 1]);
    mytitle = plotInfo.text;
    mytitle = strrep(mytitle, '_','\_');
    text(0.016, 0.97, mytitle, 'fonts', 16, 'fontw', 'bold');
end

% --- slices
for n = 1:n_slices
    ix = i_slices(n);
    % new subplot   
    %subplot(nRows, nCols, n);
    subtightplot(nRows, nCols, n, gap, marg_h, marg_w);
    hold on;    
    h_brain = image(yi, zi, squeeze(Z(ix,:,:)));
    set(h_brain, 'CData', squeeze(brain_cInds(ix,:,:))');
    h_chnls = image(yi, zi, squeeze(Z(ix,:,:)));
    set(h_chnls, 'CData', squeeze(chnls_cInds(ix,:,:))', 'AlphaData',squeeze(chnls_aData(ix,:,:))');
    
%     title(['x = ' num2str(-xi(ix))]);
    if ceil(n/nCols) == nRows && mod(n,nCols) == 1
        xlabel('x-MNI');
        ylabel('z-MNI');
    else
        set(gca, 'XTick',[], 'YTick',[]);
    end
    axis image;
    
    % click on plot to see it bigger
    set(gca, 'ButtonDownFcn','call_copy');   
    set(h_brain, 'ButtonDownFcn','call_copy');
    set(h_chnls, 'ButtonDownFcn','call_copy');
    
    % left/right orientation
    if n == 1
        txt_L = yi(1) + (yi(end)-yi(1))/100*5;        % ~ 5% offset from left  side
        txt_M = yi(1) + (yi(end)-yi(1))/100*50;       % ~ 50% offset from left  side
        txt_U = zi(end) - (zi(end)-zi(1))/100*5;      % ~ 5% offset from upper side
        txt_D = zi(1) + (zi(end)-zi(1))/100*5;        % ~ 5% offset from lower side
        text(txt_M, txt_U, 'L', 'fonts',12, 'fontw','bold', 'HorizontalAlignment','center', 'Color','w');
    end    
    if n == n_slices
        text(txt_M, txt_U, 'R', 'fonts',12, 'fontw','bold', 'HorizontalAlignment','center', 'Color','w');
    end
    % MNI coor
    text(txt_L, txt_D, ['x = ' num2str(xi(ix))], 'fonts',12, 'fontw','bold', 'HorizontalAlignment','left', 'Color','w');
end

seegv_createsavefigure(f, outDir)