function bnd_out = go_add_concentric_mesh(bnd,shift)

% Get normals of original mesh
mesh_1 = [];
mesh_1.vertices = bnd.pos;
mesh_1.faces = bnd.tri;
nrms = spm_mesh_normals(mesh_1); 

% push new vertices, then repair
mesh_2.vertices = mesh_1.vertices + shift*nrms;
mesh_2.faces = mesh_1.faces;
[mesh_2.vertices,mesh_2.faces] = meshcheckrepair(mesh_2.vertices,mesh_2.faces,'meshfix');

% combine
mesh_merge(1) = mesh_2;
mesh_merge(2) = mesh_1;
mesh_merge = spm_mesh_join(mesh_merge);

% housekeeping
bnd_out.pos = mesh_merge.vertices;
bnd_out.tri = mesh_merge.faces;
bnd_out.unit = bnd.unit;

end