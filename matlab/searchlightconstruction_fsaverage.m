
% the biggest we would want to do is probably 32-mm diameter searchlight
clear all

% setup
fsdir = [nsd_datalocation '/freesurfer'];

% define
hemis = {'lh' 'rh'};
cutoff = 10;           % all vertices within X mm of this one. the full diameter is 2*X.
dsurf = 'fsaverage3';  % fsaverage3 means every 16 mm edge distance

% do it
searchlightidxs = {};   % cell of 8 subjects x 2 hemis x L, each is row vector with indices of subject-native vertices
searchlightpeaks = {};  % cell of 8 subjects x 2 hemis, each is 1 x L with subject-native indices of where the searchlight is centered
propagatematrix = {};   % cell of 8 subjects x 2 hemis, each is V x 1 with 1-index (between 1-L) with the searchlight result to inherit

for hh=1:length(hemis)

  % load the decimated surface
  surfD = cvnreadsurface(dsurf,hemis{hh},'sphere','orig');

  % for each dsurf vertex, figure out what the corresponding fsaverage indices are (1-based)
  fsavgix = cvntransfertosubject('fsaverage',dsurf,(1:163842)',hemis{hh},'nearest','orig','orig');  % 642 x 1

    % load in inflated and sphere subject-native surfaces
    surf =  cvnreadsurface('fsaverage',hemis{hh},'inflated','orig');
    surfS = cvnreadsurface('fsaverage',hemis{hh},'sphere','orig');
    n = size(surf.vertices,1);  

    % loop over each location
    for p=1:length(fsavgix)

      % figure out which subject-native index this is
      vv = fsavgix(p);
      searchlightpeaks{1,hh}(p) = vv;  % record

      % vertex-to-vertex distance in the inflated space (3D Euclidean distance)
      dists = sqrt(sum((surf.vertices-repmat(surf.vertices(vv,:),[n 1])).^2,2));

      % record
      searchlightidxs{1,hh,p} = flatten(find(dists <= cutoff));  % all vertices within <cutoff>

    end
end
  
%% %%%%% EXPORT

for subjix=1%:8
  for hh=1:2
    searchix = [];
    for ll=1:size(searchlightidxs,3)
      searchix(ll,1:length(searchlightidxs{subjix,hh,ll})) = searchlightidxs{subjix,hh,ll};
    end
    save([nsd_datalocation('local') '/searchlight/' sprintf('searchix_radius%d_%s_fsaverage_hemi%d.mat',cutoff,dsurf,hh)],'searchix');
  end
end