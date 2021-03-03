clear all
close all

hemi = 'rh'; 
subjid = 'subj05';  subjix=5;
load('/oak/stanford/groups/kalanit/biac2/kgs/projects/Dawn/NSD/local_data/processed/subj05rh_r1corrs.mat');

% load custom map
nc = load(sprintf('/oak/stanford/groups/kalanit/biac2/kgs/projects/Dawn/NSD/local_data/freesurfer/%s/lh_split_half.mat', subjid));
lh_nc = nc.mean';
nc = load(sprintf('/oak/stanford/groups/kalanit/biac2/kgs/projects/Dawn/NSD/local_data/freesurfer/%s/rh_split_half.mat', subjid));
rh_nc = nc.mean';

i = idx'+1; %switch from python

for v = 99%:20
    temp = zeros(length(rh_nc),1);
    temp(i) = abs(matrix(v,:));

    cvnlookup(subjid,10,temp,[.001,1], jet(256), .001)
    
    pause(.2)
end

full_vi = i(v); 
surfS = cvnreadsurface(sprintf('subj%02d',subjix),hemi,'sphere','orig');
dist = temp;

for v = 1:length(surfS.vertices)
    dist(v) = norm(surfS.vertices(full_vi,:)-surfS.vertices(v,:));
end

cvnlookup(subjid,10,dist,[0,200])


full = [dist,temp];
bydist = sortrows(full);
bydist = bydist(2:end,:); %ignore self

plot(smooth(bydist(1:100000,2),100))
scatter(bydist(1:1000,1),bydist(1:100,2))
