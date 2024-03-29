%% Tessellate visual cortex into equally spaced ROIs
% this is acheived by sampling a sphere at a fixed rate and mapping this to
% the subject freesurfer sphere, where we can then assign vertices to their
% nearest sphere point

clear all
close all

%% Sample the sphere
addpath(genpath('/home/dfinzi/Desktop/S2-Sampling-Toolbox'))
sampling = 500; %300 is relatively fine 
[SP, Tri, ~,~] = ParticleSampleSphere('N', sampling); 
fv = struct('faces',Tri,'vertices',SP);
%fv = SubdivideSphericalMesh(fv,2);
SP = fv.vertices;
SP_scaled = SP;%*100;

%% Convert to specified subject and hemi
hemis = {'lh' 'rh'};

subjix=8;hh=2;
subjid='subj08';

surfS = cvnreadsurface(sprintf('subj%02d',subjix),hemis{hh},'sphere','orig');
n = size(surfS.vertices,1);

[~,iix] = max(surfS.vertices * SP_scaled', [], 2); %fast nearest neighbor equivalent

% check sizes of each parcel in num vertices
histonsph = zeros(size(SP,1),1);
for p = 1:size(SP,1)
    histonsph(p)=sum(iix==p);
end

% load stream ROI values
if hh == 1 %left
    roivals = cvnloadmgz(sprintf('/oak/stanford/groups/kalanit/biac2/kgs/projects/Dawn/NSD/data/nsddata/freesurfer/%s/label/lh.streams.mgz',subjid)); %_shrink5  % load in an existing file?
else %right
    roivals = cvnloadmgz(sprintf('/oak/stanford/groups/kalanit/biac2/kgs/projects/Dawn/NSD/data/nsddata/freesurfer/%s/label/rh.streams.mgz',subjid)); %_shrink5  % load in an existing file?
end

count = 1;
for r = 1:max(iix)  
    
    streamvals = roivals(find(iix==r));
    percent_stream = sum(streamvals~=0)/size(streamvals,1);
    
    %get rid of parcels based on overlap with stream ROIs
    if percent_stream < .5 %if less than 50% within boundaries, remove %.9
        iix(find(iix==r)) = 0;
    else %update ROI numbering
        iix(find(iix==r)) = count;
        count = count+1;
    end
end

if hh == 2 %left hemi already completed
    left = cvnloadmgz(sprintf('/oak/stanford/groups/kalanit/biac2/kgs/projects/Dawn/NSD/local_data/freesurfer/%s/lh.tessellate_500.mgz',subjid));  % load in an existing file?
    %left = cvnloadmgz(sprintf('/oak/stanford/groups/kalanit/biac2/kgs/projects/Dawn/NSD/local_data/freesurfer/%s/lh.tessellate_10vox_shrink5.mgz',subjid));  % load in an existing file?
    roivals = [left; iix];
else
    roivals = iix;  
end

% Manual plotting and saving of mgzs using cvndefinerois
roilabels=[];
cmap   = repmat(hsv(1000), 1, 1);
rng    = [0 max(roivals)];
threshs = {[] [] []};

mgznames = {'corticalsulc' 'Kastner2015', 'streams'};
crngs = {[0 28] [0 25] [0 25]};
cmaps = {jet(256) jet(256) jet(256)};
cvndefinerois;

%currently saving in local_data under ?h.tessellate_300.mgz