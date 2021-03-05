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

%% plot correlation maps and save as video
i = idx'+1; %switch from python
across = load('/oak/stanford/groups/kalanit/biac2/kgs/projects/Dawn/NSD/local_data/processed/r1r2corrs_mini.mat');

imdir = '/oak/stanford/groups/kalanit/biac2/kgs/projects/Dawn/NSD/local_data/processed/samer_images/';
for v = 1:100
    filename = fullfile(imdir, [sprintf('%03d',v) '.jpg']);
    temp = zeros(length(rh_nc),1);
    temp(i) = abs(matrix(v,:)); %max(temp(temp~=1))
    
    %across_comp= abs(across.matrix(v,:));

    [rawimg,Lookup,rgbimg] = cvnlookup(subjid,10,temp,[.001,.3], jet(256), .001);
    imwrite(rgbimg(:,1652:end,:), filename)
    
end
close all

imageNames = dir(fullfile(imdir, '*.jpg')); 
imageNames = {imageNames.name}';

writerObj = VideoWriter('r1_betacorrs.avi');
writerObj.FrameRate=3;

open(writerObj);

for ii = 1:length(imageNames)
    img = imread(fullfile(imdir, imageNames{ii})); 
    writeVideo(writerObj,img)
end
close(writerObj)

%% compute distances (on sphere as proxy)
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
scatter(bydist(1:100000,1),bydist(1:100000,2))
