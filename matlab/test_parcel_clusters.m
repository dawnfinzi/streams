clear all
close all

subjix=2;hh=2;
subjid='subj02';
sid = '02';

l = cvnloadmgz(sprintf('/oak/stanford/groups/kalanit/biac2/kgs/projects/Dawn/NSD/local_data/freesurfer/%s/lh.tessellate_1000_trim10.mgz',subjid));  % load in an existing file?
left = zeros(length(l),1);

right = cvnloadmgz(sprintf('/oak/stanford/groups/kalanit/biac2/kgs/projects/Dawn/NSD/local_data/freesurfer/%s/rh.tessellate_1000_trim10.mgz',subjid));  % load in an existing file?
%%if loading both hemis - update right roi nums to account for left
% for i = 1:length(right)
%     if right(i) ~= 0
%         right(i) = right(i) + 90;
%     end
% end
roivals = [left; right];

%% Manual plotting and saving of mgzs using cvndefinerois
load(['/oak/stanford/groups/kalanit/biac2/kgs/projects/Dawn/NSD/local_data/processed/clusters/'  subjid  '_rh_500_cluster_testing.mat']);
load(['/oak/stanford/groups/kalanit/biac2/kgs/projects/Dawn/NSD/local_data/processed/clusters/'  subjid  '_rh_500_ap.mat']);


roilabels=[];
clear col
sample_cols = viridis(single(max(full)+1)); %viridis
col = sample_cols(full+1,:);

cmap   = col; %repmat(hsv(256), 1, 1);

rng    = [1 max(right)];
threshs = {[] [] []};

mgznames = {'corticalsulc' 'Kastner2015', 'streams'};
crngs = {[0 28] [0 25] [0 25]};
cmaps = {jet(256) jet(256) jet(256)};
cvndefinerois;

%% video
load(['/oak/stanford/groups/kalanit/biac2/kgs/projects/Dawn/NSD/local_data/processed/parcel_megas/'  sid  '_1000_tt10_5000imgs_rh_mega_matrix.mat']);%matrix = squeeze(matrix);
matrix=squeeze(matrix);
extraopts = {'roiname',{'streams'},'roicolor',{'k'},'drawroinames',false, 'roiwidth', 2};

imdir = ['/oak/stanford/groups/kalanit/biac2/kgs/projects/Dawn/NSD/local_data/processed/movies/images_'  sid  '_1000_tt10/'];
if ~exist(imdir)
    mkdir(imdir)
end

for r = 1:length(matrix)
    filename = fullfile(imdir, [sprintf('%03d',r) '.jpg']);
    vec = matrix(r,:)';
    m = int16(max(vec(vec~=1))*1000);
    sample_cols = autumn(double(m)); %autumn(length(matrix));
       
    if m~=0
        [v, si] = sortrows(vec);
        idx = int16(v*1000);
        idx(idx < 1) = 1; %no zeros or negative
        idx(idx == 1000) = 1; %recolor as seed later
        full = [si, sample_cols(idx,:)];
        full = sortrows(full);
        full(r,:) = [1 0 0 0]; %seed
        cmap   = full(:,2:4);

        [rawimg,Lookup,rgbimg] = cvnlookup(subjid,13,roivals,[1,length(matrix)], double(cmap), .9,[],1,extraopts);
        imwrite(rgbimg(:,1652:end,:), filename)
    end
    close all
    
end
close all

imageNames = dir(fullfile(imdir, '*.jpg')); 
imageNames = {imageNames.name}';

cd(imdir)
writerObj = VideoWriter(['subj' sid '_1000_tt10_corrs.avi']);
writerObj.FrameRate=3;

open(writerObj);

for ii = 1:length(imageNames)
    img = imread(fullfile(imdir, imageNames{ii})); 
    writeVideo(writerObj,img)
end
close(writerObj)

