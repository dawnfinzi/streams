clear all
close all

subjix=2;hh=2;
subjid='subj02';
sid = '02';

load(['/oak/stanford/groups/kalanit/biac2/kgs/projects/Dawn/NSD/local_data/processed/parcel_megas/'  sid  '_1000_tt10_5000imgs_rh_mega_matrix.mat']);%matrix = squeeze(matrix);
matrix=squeeze(matrix);

count=1;
for r = 1:50%length(matrix)
    %filename = fullfile(imdir, [sprintf('%03d',r) '.jpg']);
    vec = matrix(r,:)';
    m = int16(max(vec(vec~=1))*1000);
    sample_cols = autumn(double(m)); %autumn(length(matrix));
       
    if m~=0
        [v, si] = sortrows(vec);
        idx = int16(v*1000);
        idx(idx < 1) = 1; %no zeros or negative
        idx(idx == 1000) = m; %1; %recolor as seed later
        full = [si, sample_cols(idx,:)];
        full = sortrows(full);
        %full(r,:) = [1 0 0 0]; %seed
        cmap   = full(:,2:4);

        [rawimg,Lookup,rgbimg] = cvnlookup(subjid,13,roivals,[1,length(matrix)], double(cmap), .9,[],1); %,extraopts);
        %imwrite(rgbimg(:,1652:end,:), filename)
        
        stack(count,:,:) = rgb2gray(double(rgbimg(:,1501:3000,:)));
        count = count + 1;
    end
    close all
    
end
close all

for s = 1:50
    edge(s,:,:) = detectedges(squeeze(stack(s,:,:)),3);
end

figure; imagesc(squeeze(mean(edge)))

%%
%%% or use r1r2betacorrs in _old/images
for n = 1:99
  
    I = imread([sprintf('%03d',n) '.jpg']);
    Ii = rgb2gray(I);
    test(n,:,:) = detectedges(double(Ii(:,1:1500)),4);
end
figure; imagesc(squeeze(mean(test)))