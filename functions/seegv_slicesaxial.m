function [ ] = seegv_slicesaxial(mni_coors, chnls_cInds,chnls_aData, ...
    brain_cInds, plotInfo, clrmap,  xi, yi, zi, Z)
%SEEGV_SLICESAXIAL Summary of this function goes here
%   Detailed explanation goes here
% --- number of slices
sliceStep = 1.5;            % in [mm]
n_maxSlices = 120;
dim = 3;
min_mni = min(mni_coors(:, dim));
max_mni = max(mni_coors(:, dim));
i_sliceStep = closestval(zi, sliceStep) - closestval(zi, 0.0);
i_slices = closestval(zi, min_mni):i_sliceStep:closestval(zi, max_mni);    
n_slices = length(i_slices);
while n_slices > n_maxSlices
    sliceStep = sliceStep + voxSize_new;        % increase slice step
    i_sliceStep = closestval(zi, sliceStep) - closestval(zi, 0.0);
    i_slices = closestval(zi, min_mni):i_sliceStep:closestval(zi, max_mni);    
    n_slices = length(i_slices);
end
% --- figure
f = figure('visible', 'on');
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
    iz = i_slices(end-n+1);
    % new subplot   
    %subplot(nRows, nCols, n);
    subtightplot(nRows, nCols, n, gap, marg_h, marg_w);
    hold on;    
    h_brain = image(xi, yi, Z(:, :, iz));
    set(h_brain, 'CData', brain_cInds(:,:,iz)');
    h_chnls = image(xi, yi, Z(:,:,iz));
    set(h_chnls, 'CData', chnls_cInds(:,:,iz)', 'AlphaData',chnls_aData(:,:,iz)');
    
%     title(['z = ' num2str(zi(iz))]);
    if ceil(n/nCols) == nRows && mod(n,nCols) == 1
        xlabel('x-MNI');
        ylabel('y-MNI');
    else
        set(gca, 'XTick',[], 'YTick',[]);
    end
    axis image;
    
    % click on plot to see it bigger
    set(gca, 'ButtonDownFcn', 'call_copy');   
    set(h_brain, 'ButtonDownFcn', 'call_copy');
    set(h_chnls, 'ButtonDownFcn', 'call_copy');
    
    % left/right orientation
    if n == 1
        txt_L = xi(1) + (xi(end)-xi(1))/100*5;        % ~ 5% offset from left  side
        txt_R = xi(end) - (xi(end)-xi(1))/100*5;      % ~ 5% offset from right side
        txt_U = yi(end) - (yi(end)-yi(1))/100*5;      % ~ 5% offset from upper side
        txt_D = yi(1) + (yi(end)-yi(1))/100*5;        % ~ 5% offset from lower side
        text(txt_L, txt_U, 'R', 'fonts',12, 'fontw','bold', 'HorizontalAlignment','center', 'Color','w');
        text(txt_R, txt_U, 'L', 'fonts',12, 'fontw','bold', 'HorizontalAlignment','center', 'Color','w');
    end
    % MNI coor
    text(txt_L, txt_D,['z = ' num2str(zi(iz))], 'fonts',12, 'fontw','bold', 'HorizontalAlignment','left', 'Color','w');
end

seegv_createsavefigure(f, outDir)
