clear all
close all

subjix=6;hh=2;
subjid='subj06';

l = cvnloadmgz(sprintf('/oak/stanford/groups/kalanit/biac2/kgs/projects/Dawn/NSD/local_data/freesurfer/%s/lh.tessellate_500.mgz',subjid));  % load in an existing file?
left = zeros(length(l),1);

right = cvnloadmgz(sprintf('/oak/stanford/groups/kalanit/biac2/kgs/projects/Dawn/NSD/local_data/freesurfer/%s/rh.tessellate_500.mgz',subjid));  % load in an existing file?
%%if loading both hemis - update right roi nums to account for left
% for i = 1:length(right)
%     if right(i) ~= 0
%         right(i) = right(i) + 90;
%     end
% end
roivals = [left; right];

%% Manual plotting and saving of mgzs using cvndefinerois
roilabels=[];
clear col
sample_cols = viridis(15);
col = sample_cols(full+1,:);

% sample_cols = jet(5);
% for r = 1:max(roivals);
%     if split1(r) == 0
%         col(r,:) = [1, 1, 1];
%     elseif split1(r) == 1
%         col(r,:) = sample_cols(1,:);
%     end
% end
% for r = 1:max(roivals);
%     if split4(r) == 0
%         col(r,:) = [1, 1, 1];
%     elseif split4(r) == 1
%         col(r,:) = sample_cols(1,:);
%     elseif split4(r) == 2
%         col(r,:) = sample_cols(2,:);
%     elseif split4(r) == 3
%         col(r,:) = sample_cols(3,:);
%     elseif split4(r) == 4
%         col(r,:) = sample_cols(4,:);
%     elseif split4(r) == 5
%         col(r,:) = sample_cols(5,:);
%     elseif split4(r) == 6
%         col(r,:) = sample_cols(6,:);
%     elseif split4(r) == 7
%         col(r,:) = sample_cols(7,:);
%     elseif split4(r) == 8
%         col(r,:) = sample_cols(8,:);
%     elseif split4(r) == 9
%         col(r,:) = sample_cols(9,:);
%     elseif split4(r) == 10
%         col(r,:) = sample_cols(10,:);
%     end
% end
% sample_cols = hsv(10);
% for r = 1:max(roivals);
%     if ap(r) == 0
%         col(r,:) = sample_cols(8,:);
%     elseif ap(r) == 1
%         col(r,:) = sample_cols(1,:);
%     elseif ap(r) == 2
%         col(r,:) = sample_cols(2,:);
%     elseif ap(r) == 3
%         col(r,:) = sample_cols(3,:);
%     elseif ap(r) == 4
%         col(r,:) = sample_cols(4,:);
%     elseif ap(r) == 5
%         col(r,:) = sample_cols(5,:);
%     elseif ap(r) == 6
%         col(r,:) = sample_cols(6,:);
%     elseif ap(r) == 7
%         col(r,:) = sample_cols(7,:);
%     elseif ap(r) == 8
%         col(r,:) = sample_cols(8,:);
%     elseif ap(r) == 9
%         col(r,:) = sample_cols(9,:);
%     end
% end

cmap   = col; %repmat(hsv(256), 1, 1);

rng    = [1 max(right)];
threshs = {[] [] []};

mgznames = {'corticalsulc' 'Kastner2015', 'streams'};
crngs = {[0 28] [0 25] [0 25]};
cmaps = {jet(256) jet(256) jet(256)};
cvndefinerois;

%% video
matrix = squeeze(matrix);
extraopts = {'roiname',{'streams'},'roicolor',{'k'},'drawroinames',false, 'roiwidth', 2};

imdir = '/oak/stanford/groups/kalanit/biac2/kgs/projects/Dawn/NSD/local_data/processed/images_06_500/';
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
    
end
close all

imageNames = dir(fullfile(imdir, '*.jpg')); 
imageNames = {imageNames.name}';

writerObj = VideoWriter('subj06_500_scale_corrs.avi');
writerObj.FrameRate=3;

open(writerObj);

for ii = 1:length(imageNames)
    img = imread(fullfile(imdir, imageNames{ii})); 
    writeVideo(writerObj,img)
end
close(writerObj)

