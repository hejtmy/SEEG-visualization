function [ brain ] = seegv_interpolate( brain, interpolatedSize )
%SEEGV_INTERPOLATE Summary of this function goes here
%   Detailed explanation goes here
% brain: created by getbraindata.m - spm brain created by spmvol and smp_read_vols
    %% MNI axis (in mm)
    voxSize_old = [NaN, NaN, NaN];
    tmp = reshape(brain.xyz(1, :), size(brain.vol));
    x = -squeeze(tmp(:, 1, 1));
    assert(x(1) < x(end));          % ascending order
    voxSize_old(1) = x(2) - x(1);

    tmp = reshape(brain.xyz(2, :), size(brain.vol));
    y = squeeze(tmp(1, :, 1))';
    assert(y(1) < y(end));
    voxSize_old(2) = y(2) - y(1);

    tmp = reshape(brain.xyz(3, :), size(brain.vol));
    z = squeeze(tmp(1, 1, :));
    assert(z(1) < z(end));
    voxSize_old(3) = z(2) - z(1);

    %% interpolation (if needed)
    voxSize_new = interpolatedSize;     % in [mm]
    if all(voxSize_new == voxSize_old)
        brainVolume = brain.vol;
        xi = x';
        yi = y';
        zi = z';
    else
        [X, Y, Z] = meshgrid(x, y, z);
        xi = x(1):voxSize_new:x(end);
        yi = y(1):voxSize_new:y(end);
        zi = z(1):voxSize_new:z(end);
        [XI, YI, ZI] = meshgrid(xi, yi, zi);
        V = permute(brain.vol, [2, 1, 3]);                        % permutation needed by interp3 (for some reason)
        brainVolume = interp3(X, Y, Z, V, XI, YI, ZI, 'spline');
        brainVolume = permute(brainVolume, [2, 1, 3]);                                   % re-arrange back to [x,y,z] dimensions
    end

    %% return interpolated volume and axis
    assert(size(brainVolume, 1) == size(xi, 2));
    assert(size(brainVolume, 2) == size(yi, 2));
    assert(size(brainVolume, 3) == size(zi, 2));
    brain.VI = brainVolume;
    brain.xi = xi;
    brain.yi = yi;
    brain.zi = zi;
    brain.voxSize_new = voxSize_new;

end
