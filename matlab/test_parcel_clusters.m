clear all
close all

subjix=5;hh=2;
subjid='subj05';

l = cvnloadmgz(sprintf('/oak/stanford/groups/kalanit/biac2/kgs/projects/Dawn/NSD/local_data/freesurfer/%s/lh.tessellate_300.mgz',subjid));  % load in an existing file?
left = zeros(length(l),1);

right = cvnloadmgz(sprintf('/oak/stanford/groups/kalanit/biac2/kgs/projects/Dawn/NSD/local_data/freesurfer/%s/rh.tessellate_300.mgz',subjid));  % load in an existing file?
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
sample_cols = hsv(20);

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

cmap   = col;
%cmap   = repmat(hsv(256), 1, 1);

rng    = [0 max(right)+1];
threshs = {[] [] []};

mgznames = {'corticalsulc' 'Kastner2015', 'streams'};
crngs = {[0 28] [0 25] [0 25]};
cmaps = {jet(256) jet(256) jet(256)};
cvndefinerois;