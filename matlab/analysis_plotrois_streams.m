%% when done with all subjects, write out inflated views

% these images are our final inspections and quality assessments
clear all
outputdir = '/oak/stanford/groups/kalanit/biac2/kgs/projects/Dawn/NSD/results/ROIs';
mkdirquiet(outputdir);

%%
% let's look at floc and stream ROIs at the same time
for subjix=1:8
  extraopts = {'roiname',{'floc-faces' 'floc-words' 'streams'},'roicolor',{'r' 'b' 'k'},'drawroinames',true, 'roiwidth', 2};
  cvnlookup(sprintf('subj%02d',subjix),13,[],[],[],100,[],1,extraopts);
  imwrite(rgbimg,sprintf('%s/subj%02d_floc-streams-names.png',outputdir,subjix));
  
  extraopts = {'roiname',{'floc-faces' 'floc-words' 'floc-places' 'floc-bodies' 'streams'},'roicolor',{'r' 'b' 'g' 'y' 'k'},'drawroinames',false, 'roiwidth', 1.5};
  cvnlookup(sprintf('subj%02d',subjix),13,[],[],[],100,[],1,extraopts);
  imwrite(rgbimg,sprintf('%s/subj%02d_all-floc-streams.png',outputdir,subjix));
end

%% now look at eccentricity and phase maps with stream ROIs
for subjix=1:8
   % load ecc map
   file0 = [nsd_datalocation sprintf('/freesurfer/subj%02d/label/*.prfeccentricity.mgz',subjix)];
   file0 = cvnloadmgz(file0); 
   
   ve_mask = [nsd_datalocation sprintf('/freesurfer/subj%02d/label/*.prfR2.mgz',subjix)];
   ve_mask = cvnloadmgz(ve_mask); 
   
   file0(ve_mask<10) = 0;
    
  extraopts = {'roiname',{'streams'},'roicolor',{'k'},'drawroinames',true, 'roiwidth', 2};
  cvnlookup(sprintf('subj%02d',subjix),13,file0,[0.01 7.5],colormap(flipud(jet)),0.1,[],1,extraopts);
 % imwrite(rgbimg,sprintf('%s/subj%02d_ecc-streams.png',outputdir,subjix));
end

%phase
for subjix=1:8
   % load phase map
   file0 = [nsd_datalocation sprintf('/freesurfer/subj%02d/label/*.prfangle.mgz',subjix)];
   file0 = cvnloadmgz(file0); 
   
   ve_mask = [nsd_datalocation sprintf('/freesurfer/subj%02d/label/*.prfR2.mgz',subjix)];
   ve_mask = cvnloadmgz(ve_mask); 
   
   file0(ve_mask<10) = 0;
   
   cmap_phase = load('WedgeMapRight_pRF.mat');
   cmr = cmap_phase.modeInformation.cmap;
   cmap_phase = load('WedgeMapLeft_pRF.mat');
   cml = cmap_phase.modeInformation.cmap;
   cm = [cmr(75:118,:); cmr(75:118,:)];
    
  extraopts = {'roiname',{'streams'},'roicolor',{'k'},'drawroinames',true, 'roiwidth', 2};
  cvnlookup(sprintf('subj%02d',subjix),13,file0,[0.001 361],cm,0.001,[],1,extraopts);
  imwrite(rgbimg,sprintf('%s/subj%02d_phase-streams.png',outputdir,subjix));
end

%% visualize noise ceiling maps too
for subjix=1:8
    lh_a1 = load_mgh([sprintf('%s/ppdata/subj%02d/nativesurface/betas_fithrf_GLMdenoise_RR/',nsd_datalocation('betas'),subjix) 'lh.nc_3trials.mgh']);  
    rh_a1 = load_mgh([sprintf('%s/ppdata/subj%02d/nativesurface/betas_fithrf_GLMdenoise_RR/',nsd_datalocation('betas'),subjix) 'rh.nc_3trials.mgh']);  
    
    file0 = [lh_a1; rh_a1];
    
    extraopts = {'roiname',{'streams'},'roicolor',{'k'},'drawroinames',true, 'roiwidth', 2};
    cvnlookup(sprintf('subj%02d',subjix),13,file0,[0 75],jet,20,[],1,extraopts);
    imwrite(rgbimg,sprintf('%s/subj%02d_noiseceiling-streams.png',outputdir,subjix));
end

%% same for split half maps _FIX!!
for subjix=1:8

    lh_a1 = load_mgh([nsd_datalocation('local'), '/freesurfer/subj0' , num2str(subjix), '/lh.avg_split_half.mgh']);  
    rh_a1 = load_mgh([nsd_datalocation('local'), '/freesurfer/subj0' , num2str(subjix), '/rh.avg_split_half.mgh']);  
    
    file0 = [lh_a1; rh_a1];
    
    extraopts = {'roiname',{'streams'},'roicolor',{'k'},'drawroinames',true, 'roiwidth', 2};
    cvnlookup(sprintf('subj%02d',subjix),13,file0,[0 .5],jet,0.05,[],1,extraopts);
    imwrite(rgbimg,sprintf('%s/subj%02d_split-half-streams.png',outputdir,subjix));
end
