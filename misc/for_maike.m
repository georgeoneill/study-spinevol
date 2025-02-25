clearvars
close all
clc

proj_init

load('C:\Users\George O''Neill\Downloads\geometries_19seg_bem.mat')

bv = spm_mesh_bounding_volume(mesh_wm);

tmp = [];
tmp.p = mesh_wm.vertices;
tmp.e = mesh_wm.faces;

tmp = hbf_CorrectTriangleOrientation(tmp);


% my function uses the fieldtrip naming convention
mesh2.pos = tmp.p;
mesh2.tri = tmp.e;
mesh2.unit = 'm';


cent = [0, -0.027, -0.110];

% Now we hand over to iso2mesh to generate the tetrahedral mesh, I think we
% have told it no tetrahedron to be no larger than 10 ml
[node,elem] = surf2mesh(mesh2.pos, mesh2.tri, bv(1,:), bv(2,:),...
    1, 5*1e-7, cent,[],[],'tetgen1.5');
fem_plot_mesh(node,elem,3,'b',[-0.02 0])
axis equal

% make a mesh 2mm inwards and add to original
mesh2 = go_add_concentric_mesh(mesh2,0.002);
% tessalate
[node,elem] = surf2mesh(mesh2.pos, mesh2.tri, bv(1,:), bv(2,:),...
    1, 5*1e-7, cent,[],[],'tetgen1.5');
elem(:,end) = 1;
fem_plot_mesh(node,elem,3,'b',[-0.02 0])
axis equal