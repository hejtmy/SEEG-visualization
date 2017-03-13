% requires only plotInfo about file and interpolation size. Jsut pass it as parameters
function brain = getbraindata(filename)
% loads 3D brain MRI or CT = brain
% interpolates
% finds selected channels to plot
% MUST inlcude SPM package in Matlab path

% (c) Jiri, Jan17

%% is SPM12 installed ?
spm_dir = what('spm12');
if isempty(spm_dir) 
    error('SPM12 toolbox required. Install SPM12 and add it to pathdef.m !'); 
end

%% load normalized brain MRI (wT1.nii)
assert(exist(filename, 'file') == 2);
brain.hdr = spm_vol(filename);
[brain.vol, brain.xyz] = spm_read_vols(brain.hdr);
