function tet = go_paperbnds2mesh_denserspine(bnd,nshells)
% A version of the FEM meshing code which incorporates the fact the WM mesh
% has been turned into a set of concentric shells
%% Turn the boundaries into one tesselated mesh 
%
% 1) Define seed locations which are inside a given boundary 
%
% 2) Take the boundaries from a mesh and turn it into a tetrahedral mesh for
% finite element analyses

if nargin == 1
    error(['please specify how many shells ',...
        'the spinal cord has been subvided into!'])
end

% First merge all boundaries into one mesh...
bemMerge = {};
for ii = 1:numel(bnd)
    bemMerge = cat(2, bemMerge, bnd(ii).pos, bnd(ii).tri);
end
[newnode, newelem] = mergemesh(bemMerge{:});
tmp = [];
tmp.vertices = newnode;
tmp.faces = newelem(:,1:3);

% ...then break it up again...
klust = spm_mesh_split(tmp);
% ...and merge it again!
bemMerge = {};
for ii = 1:numel(klust)
    bemMerge = cat(2, bemMerge, klust(ii).vertices, klust(ii).faces);
end
[newnode, newelem] = mergemesh(bemMerge{:});

% For the non-spine parts of the body, just find a location inside to
% generate seed locations (as these are convex hulls the mean vertex is 
% likely to work) We split earlier as we wanted to treat each lung and
% ventricle seperately
for ii = (nshells+2):numel(klust)
    tmp = klust(ii).vertices;
    cent(ii,:) = mean(tmp);
end

% Working out where to find a seed for the spine compartments, as these are
% narrow and wriggly we should be a little more careful finding points
% which are inside. I think I've defined a plane roughly halfway up the
% spine and we just slowly move outwards
fprintf('determing point inside innermost mesh...')
bv = spm_mesh_bounding_volume(klust(1));
[~, dimmax] = max(abs(diff(bv)));
rng{1} = bv(1,1):0.001:bv(2,1);
rng{2} = bv(1,2):0.001:bv(2,2);
rng{3} = bv(1,3):0.001:bv(2,3);
rng{dimmax} = 0.5*(bv(1,dimmax) + bv(2,dimmax));
[xx, yy, zz] = ndgrid(rng{1},rng{2},rng{3});
candidates = [xx(:) yy(:) zz(:)];
inside = zeros(length(candidates),1);
for ii = 1:length(candidates)
    inside(ii) = tt_is_inside(candidates(ii,:),...
        klust(1).vertices,klust(1).faces);
end
cent(1,:) = mean(candidates(find(inside),:));
assert(tt_is_inside(cent(1,:),klust(1).vertices,klust(1).faces),...
    'random seed source not inside spinal column!');
fprintf(' COMPLETE\n')

% shift it forward to work out when we have left one mesh at a time up to
% the bone
fprintf('determing seeds for other meshes...\n')
shifts = (1:30)' * [0.0005 0 0];
points = cent(1,:) + shifts;
inside = zeros(nshells+1,length(points));
for ii = 1:nshells+1
    for jj = 1:length(points)
        inside(ii,jj) = tt_is_inside(points(jj,:),...
            klust(ii).vertices,klust(ii).faces);
    end
end
inshells = sum(inside);
% get the seeds to give to the meshes
for ii = 1:nshells+1
    id = find(inshells == 2 + nshells - ii);
    cent(ii,:) = mean(points(id,:));
end
fprintf(' COMPLETE\n')

% Now we hand over to iso2mesh to generate the tetrahedral mesh, I think we
% have told it no tetrahedron to be no larger than 10 ml
[node,elem] = surf2mesh(newnode, newelem, min(newnode), max(newnode),...
    1, 5*1e-7, cent,[],[],'tetgen1.5');

% There are n compartments IDs, but we only have 5 tissue types, so this
% fixes it
% Turn first nshells ID into 1 
id = find(elem(:,5) <= nshells);
elem(id,5) = 1;
% now do the bone which should be equal to nshells + 1;
id = find(elem(:,5) == nshells+1);
elem(id,5) = 2;
% next compariment was split into two so looking for + 2 and +3
id = find(elem(:,5) == nshells+2 | elem(:,5) == nshells+3);
elem(id,5) = 3;
% next compariment was split into two so looking for + 4 and +5
id = find(elem(:,5) == nshells+4 | elem(:,5) == nshells+5);
elem(id,5) = 4;
% then the last one should be nshells + 6
id = find(elem(:,5) == nshells+6);
elem(id,5) = 5;

% Some general housekeeping to ensure the mesh itegrity
[no,el] = removeisolatednode(node,elem(:,1:4));
newelem = meshreorient(no, el(:,1:4));
elem = [newelem elem(:,5)];
node = no;

% Lets plot the mesh and take a slice through it
fem_plot_mesh(node,elem,3,'<',0);
axis equal
axis off
set(gcf,'color','w')

% Store as a structure for later
tet = [];
tet.pos = node;
tet.tet = elem(:,1:4);
tet.tissue = elem(:,5);
tet.unit = 'm';