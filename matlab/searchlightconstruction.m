%%%% This is sort of "in progress". Created for Tom.

% the biggest we would want to do is probably 32-mm diameter searchlight

% setup
fsdir = [nsd_datalocation '/freesurfer'];

% define
hemis = {'lh' 'rh'};
cutoff = 12;           % all vertices within X mm of this one. the full diameter is 2*X.
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

  % loop over subjects
  for subjix=1:8

    % load in inflated and sphere subject-native surfaces
    surf =  cvnreadsurface(sprintf('subj%02d',subjix),hemis{hh},'inflated','orig');
    surfS = cvnreadsurface(sprintf('subj%02d',subjix),hemis{hh},'sphere','orig');
    n = size(surf.vertices,1);

    % for each fsaverage vertex, inherit the nearest subject-native vertex (1-based)
    subjnative = nsd_mapdata(subjix,sprintf('%s.white',hemis{hh}),'fsaverage',(1:n)','nearest');  % 163842 x 1

    % loop over each location
    for p=1:length(fsavgix)

      % figure out which subject-native index this is
      vv = subjnative(fsavgix(p));
      searchlightpeaks{subjix,hh}(p) = vv;  % record

      % vertex-to-vertex distance in the inflated space (3D Euclidean distance)
      dists = sqrt(sum((surf.vertices-repmat(surf.vertices(vv,:),[n 1])).^2,2));

      % record
      searchlightidxs{subjix,hh,p} = flatten(find(dists <= cutoff));  % all vertices within <cutoff>

    end

%     % figure out the propagation of values
%     ix = searchlightpeaks{subjix,hh};  % where do we have searchlights centered?
%     f = griddata(surfS.vertices(ix,1),surfS.vertices(ix,2),surfS.vertices(ix,3),(1:length(ix))', ...
%                  surfS.vertices(:,1),surfS.vertices(:,2),surfS.vertices(:,3),'nearest');  % V x 1 with index of nearest searchlight
%     propagatematrix{subjix,hh} = f;

  end
  
end

% % save
% save([nsd_datalocation('local') '/searchlight/' sprintf('searchix_radius%d_%s.mat',cutoff,dsurf)],'cutoff','dsurf','searchlightidxs','searchlightpeaks','propagatematrix');
% 
% % propagatematrix{subjix,hemi} has V x 1 where elements are the 1-indexed searchlight from which values should be propagated

%% %%%%% EXPORT

for subjix=1:8
  for hh=1:2
    searchix = [];
    for ll=1:size(searchlightidxs,3)
      searchix(ll,1:length(searchlightidxs{subjix,hh,ll})) = searchlightidxs{subjix,hh,ll};
    end
    save([nsd_datalocation('local') '/searchlight/' sprintf('searchix_radius%d_%s_subj%02d_hemi%d.mat',cutoff,dsurf,subjix,hh)],'searchix');
  end
end

% searchix_radius16_fsaverage3_subj01_hemi1.mat has:
%   <searchix> is M x N matrix. M refers to different searchlights.
%     N is the maximum number of vertices involved in any searchlight.
%     Each row contains the 1-based indices of the surface vertices that are
%     involved in that searchlight. Note that the two hemispheres are treated
%     separately, and the 1-based indexing restarts for each hemisphere.
%     Because different searchlights have different numbers of vertices,
%     zeros (0) are used to indicate "no data here".

%% %%%%%% NOTES

% - doing on fsaverage3 makes the visualization very discrete. should we interpolate or somehow smooth?

%%%%% JUNK

% ok = [];
% for ll=1:642
%   ok(ll) = mean(a1(searchlightidxs{2,1,ll}));
% end
% 
% fulldata = ok(propagatematrix{2,1});