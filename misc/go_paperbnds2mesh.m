function tet = go_bnds2mesh_maike(bnd,centroids,tissue);
% A version of the FEM meshing code which incorporates the fact the WM mesh
% has been turned into a set of concentric ssh
%% Turn the boundaries into one tesselated mesh 
%
% 1) Define seed locations which are inside a given boundary 
%
% 2) Take the boundaries from a mesh and turn it into a tetrahedral mesh for
% finite element analyses

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
for ii = 3:numel(klust)
    tmp = klust(ii).vertices;
    cent(ii,:) = mean(tmp);
end

% Working out where to find a seed for the spine compartments, as these are
% narrow and wriggly we should be a little more careful finding points
% which are inside. I think I've defined a plane roughly halfway up the
% spine and we just slowly move outwards
box_min = min(klust(1).vertices);
box_max = max(klust(1).vertices);
[~, dimmax] = max(abs(box_max-box_min));
rng{1} = box_min(1):0.005:box_max(1);
rng{2} = box_min(2):0.005:box_max(2);
rng{3} = box_min(3):0.005:box_max(3);
rng{dimmax} = 0.5*(box_max(dimmax) + box_min(dimmax));
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

% shift it forward untill not inside 1 but inside 2
tmp = cent(1,:);
isin1 = true;
isin2 = true;
shift1 = 0;
while isin1 && isin2
    shift1= shift1 + [0.001 0 0];
    isin1 = tt_is_inside(tmp + shift1,klust(1).vertices,klust(1).faces);
    isin2 = tt_is_inside(tmp + shift1,klust(2).vertices,klust(2).faces);
end
shift2 = shift1;
while isin2
    shift2= shift2 + [0.001 0 0];
    isin2 = tt_is_inside(tmp + shift2,klust(2).vertices,klust(2).faces);
end
shift_final = 0.5*(shift1 + shift2);
cent(2,:) = cent(1,:)+shift_final;

% Now we hand over to iso2mesh to generate the tetrahedral mesh, I think we
% have told it no tetrahedron to be no larger than 10 ml
[node,elem] = surf2mesh(newnode, newelem, min(newnode), max(newnode),...
    1, 5*1e-7, cent,[],[],'tetgen1.5');

% There are 7 compartments IDs, but we only have 5 tissue types, so this
% fixes it
id = elem(:,5) + 10;
id(id==11) = 1;
id(id==12) = 2;
id(id==13 | id == 14) = 3;
id(id==15 | id == 16) = 4;
id(id==17) = 5;

elem(:,end) = id;

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