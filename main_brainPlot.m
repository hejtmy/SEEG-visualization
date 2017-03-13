%% master script for plotting SEEG values in slices and in 3D model
% this script is meant as a tutorial showing a with 3 simple examples
% install: 
%   - add this script to matlab paths (use: Set Paths)
%   - SPM12 (or SPM8) must be installed (add paths to SPM)
% data requirements: 
%   - subject specific MRI (T1.nii) normalized to MNI space
%   - MNI coordinates of individual channels
%   - values of individual channels
% does the following:
%   (A) loads universal brain template (colin 27) and plots data in slices
%   (B) 3D brain model (semi-transparent) with colored channel values
%   (C) loads subject-specific brain data and plots them in slices

% (c) Jiri Hammer, Jan 2017, bug reports: jirihammer@gmail.com

%% set path to this script
% BAD practice 
% - if developing for this package, should be loaded 
% - if used overall, it shoudl be in the path
% - otherwise, you load it in appropriate screen, no other stuff needed

dir_toolbox = pwd;
%dir_toolbox = what('SEEG-visualization');     % folder 'visualization_SEEG' MUST be in pathdef.m
% if isempty(dir_toolbox) 
%     error('Path to toolbox required. Add it to pathdef.m !'); 
% end
% addpath(genpath(dir_toolbox.path));
% dir_curr = pwd;
% cd(dir_toolbox.path);

%% ######################## USER INTERFACE ##################################
% --- interface structure 'plotInfo'
plotInfo = seegv_skeletoncfg;

% --- load channels MNI coors (variable 'data_channels' in 'channelsInfo.mat')
load('channelsInfo.mat', 'data_channels');
assert(exist('data_channels','var') == 1);
plotInfo.chnls = data_channels;             % !!! set here your MNI coordinates as a structure array with the same fields as in this example!

% --- channels values to plot, format = [channels x time]
vals = randn(size(data_channels,2),2);      % !!! set here your channel valus, format = [channels x time], in this example: 2 time points of channel values

% --- channels color scale
if isempty(plotInfo.colorScale)
    clims = [prctile(vals(:),5), prctile(vals(:),95)]; clims = max(abs(clims)); 
    plotInfo.chnl_clims = [-clims, clims];
else
    plotInfo.chnl_clims = plotInfo.colorScale;
end



%% ################## (A) COLIN27 BRAIN SLICES ##########################
% colin27-specific: get brain template for slices
plotInfo.plottingStyle = 'slices';
plotInfo.MRI_file = [plotInfo.MRI_fileDir filesep 'wT1_subject.nii'];       % subject specific brain normalized to MNI space
plotInfo.brain = getBrainData(plotInfo);

% colin27-specific: plot time-series of brain slices (1 figure / time point)
plotInfo.outDir = [plotInfo.outputDir filesep 'slices_normalizedBrain'];
for t = 1:size(vals,2)          % go thru all time points
    plotInfo.figName = [plotInfo.figureNamePrefix 'colin_time' num2str(t)];
    plot_brainSlices(vals(:,t), plotInfo);
end    



%% ################## (B) 3D BRAIN MODEL ##########################
% colin27-specific (but works on other brain's too! Requires segmented grey matter, normalized to MNI space)
% get brain template for 3D model
plotInfo.plottingStyle = '3D_model';
plotInfo.size_interpolate = 1;          % consider larger voxels for faster & easier manipulation
plotInfo.MRI_file = [plotInfo.MRI_fileDir filesep  'wc1T1_colin27.nii'];       % gray matter (segmented from colin27 by SPM12)
plotInfo.brain = getBrainData(plotInfo);

% colin27-specific: plot brain: 3D model 
plotInfo.outDir = [plotInfo.outputDir filesep '3D_model'];
for t = 1:size(vals,2)          % go thru all time points
    plotInfo.figName = [plotInfo.figureNamePrefix 'colin_time' num2str(t)];
    plot_brain3D(vals(:,t), plotInfo);
end



%% ################## (C) SUBJECT-SPECIFIC SLICES ##########################
% subject-specific: get brain info for slices
plotInfo.plottingStyle = 'slices';
plotInfo.MRI_file = [plotInfo.MRI_fileDir filesep 'wT1_subject.nii'];       % subject specific brain normalized to MNI space
plotInfo.brain = getBrainData(plotInfo);

% subject-specific: plot time-series of brain slices (1 figure / time point)
plotInfo.outDir = [plotInfo.outputDir filesep 'slices_subjectBrain'];
for t = 1:size(vals,2)          % go thru all time points
    plotInfo.figName = [plotInfo.figureNamePrefix 'subject_time' num2str(t)];
    plot_brainSlices(vals(:,t), plotInfo);
end    


%% change back to previous directory
display('Done. Thanks for trying this tutorial out!  :o) ');
cd(dir_curr);
