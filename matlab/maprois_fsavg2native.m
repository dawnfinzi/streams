%% Map ROIs drawn on fsaverage to native surface spaces via nearest-neighbor

clear all

hemis = {'lh', 'rh'}; 

fsbase = [nsd_datalocation '/freesurfer/subj%02d'];
fsavg = [nsd_datalocation '/freesurfer/fsaverage'];

for subjix = 1:8 %eight subjects
    for p = 1:length(hemis)
        
        sourcedata = sprintf([fsavg '/label/%s.streams.mgz'], hemis{p});
        
        fsdir = sprintf(fsbase, subjix);
        outfile = sprintf([fsdir '/label/%s.streams.mgz'], hemis{p});
    
        nsd_mapdata(subjix,'fsaverage',sprintf('%s.white',hemis{p}),sourcedata,[],0,outfile,[],fsbase);
    end
    
end
