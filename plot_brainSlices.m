function plot_brainSlices(vals, outDir, plotInfo)
% plots topologically colorcoded values 'vals' of channels onto brain slices
% (c) Jiri, Jan17

%% color maps
clrmap.brain = gray(128);                           % colormap for brain (T1, T2, CT, ...)

if isfield(plotInfo, 'colorMap')
    clrmap.chnls = plotInfo.colorMap;
else
    clrmap.chnls = jet(128);                            % default colormap for values: jet
    % clrmap.chnls = getColorMap('bwr', 128);             % colormap for values: blue - white - red
    % clrmap.chnls = getColorMap('bcwwmr', 128);          % colormap for values: blue - cyan - white - magenta - red
end
clrmap.fig = cat(1, clrmap.brain, clrmap.chnls);    % colormap of the figure
alphaVal = 1.0;                                     % transparency of the colored values (1 = opaque)

%% pass info from loaded brain (see getBrainData.m)
VI = plotInfo.brain.VI;      % interpolated volume
xi = plotInfo.brain.xi;      % interpolated x-axis, in [mm] of MNI coors
yi = plotInfo.brain.yi;      % interpolated y-axis, in [mm] of MNI coors
zi = plotInfo.brain.zi;      % interpolated z-axis, in [mm] of MNI coors
assert(size(plotInfo.chnls, 2) == size(vals, 1));
voxSize_new = plotInfo.brain.voxSize_new;

%% area for colored channel values: enlarged voxel
enlargedVoxel.size_mm = plotInfo.size_coloredCube;                             % size of the enlarged "voxel", in [mm]
enlargedVoxel.size_ind = abs(closestval(xi,enlargedVoxel.size_mm/2) - closestval(xi,0));    % half-size of the enlarged "voxel", in [indices]
if enlargedVoxel.size_ind == 0, enlargedVoxel.size_ind = 1; end
enlargedVoxel.side = -enlargedVoxel.size_ind:enlargedVoxel.size_ind;    % side of the cube, in [indices], w.r.t. its center

%% insert channel values at MNI coors to 3D arrays: color and transparency
chnls_cData = zeros(size(VI));                          % color data values
chnls_aData = zeros(size(VI));                          % fully transparent
mni_coors = [];
for ch = 1:size(plotInfo.chnls, 2)
    mni_coors = cat(1, mni_coors, [-plotInfo.chnls(ch).MNI_x, plotInfo.chnls(ch).MNI_y, plotInfo.chnls(ch).MNI_z]);
    [ix, iy, iz] = mni2vox(-plotInfo.chnls(ch).MNI_x, plotInfo.chnls(ch).MNI_y, plotInfo.chnls(ch).MNI_z, xi, yi, zi); % index of MNI coor
    
    i_sel_x = ix + enlargedVoxel.side;                  % x-indices of enlarged voxel
    i_sel_x(i_sel_x < 1) = [];                          % indices out of range
    i_sel_x(i_sel_x > length(xi)) = [];
    
    i_sel_y = iy + enlargedVoxel.side;                  % y-indices of enlarged voxel
    i_sel_y(i_sel_y < 1) = [];                          % indices out of range
    i_sel_y(i_sel_y > length(yi)) = [];    
    
    i_sel_z = iz + enlargedVoxel.side;                  % z-indices of enlarged voxel
    i_sel_z(i_sel_z < 1) = [];                          % indices out of range
    i_sel_z(i_sel_z > length(zi)) = [];    
    
    chnls_cData(i_sel_x, i_sel_y, i_sel_z) = vals(ch);    % insert the chnl value in a larger (cubic) voxel
    chnls_aData(i_sel_x, i_sel_y, i_sel_z) = alphaVal;    % transparency of the chnl value in a larger (cubic) voxel
end

%% map values to colormaps
% 3D brain data
%brain_cInds = cVals2cInds(VI, [min(VI(:)),max(VI(:))], [1,size(clrmap.brain,1)]);
brain_cInds = cVals2cInds(VI, [prctile(VI(:), 1), prctile(VI(:), 99)], [1, size(clrmap.brain, 1)]);

% channel values
if ~isfield(plotInfo, 'chnl_clims')
    clims = [min(chnls_cData(:)), max(chnls_cData(:))];
else
    clims = plotInfo.chnl_clims;
end
chnls_cInds = cVals2cInds(chnls_cData, [clims(1),clims(2)], size(clrmap.brain,1)+[1,size(clrmap.chnls,1)]);
Z = ones(size(VI));     % dummy variable for image

%% ------------------------- AXIAL SLICES ---------------------------------
if ismember('axial', plotInfo.slicePlanes)
    f = seegv_slicesaxial(mni_coors, chnls_cInds, chnls_aData, brain_cInds, plotInfo, clrmap,  xi, yi, zi, Z);
    seegv_createsavefigure(f, clims, clrmap, 'slices')
end

%% ------------------------- SAGITTAL SLICES ---------------------------------
if ismember('sagittal', plotInfo.slicePlanes)
    seegv_slicessagital(mni_coors, chnls_cInds);
end
%% ------------------------- CORONAL SLICES ---------------------------------
if ismember('coronal', plotInfo.slicePlanes)
    seegv_slicescoronal(mni_coors, chnls_cInds);
end