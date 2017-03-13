function [ fv ] = seegv_isosurface( brainVolume )
%SEEGV_ISOSURFACE Summary of this function goes here
%   Detailed explanation goes here
 %% isosurface for 3D model
    V = linTransform(brainVolume, [min(brainVolume(:)), max(brainVolume(:))], [0, 1]);
    separationThreshold = 0.5;                  % (value 0.5 separates gray matter from the dark background). Other values 0 - 1 may work also fine
    fv = isosurface(V, separationThreshold);    % surface, vertex 

end

