function fem_write_gmsh(file,mesh)

% Reads in mesh structure and outputs a .msh file for use with gmsh
% applications, with a particular empahsis to get working with duneuro.

% Based on source code at:
% https://people.sc.fsu.edu/~jburkardt/m_src/gmsh_io/gmsh_io.html

% First thing to check is whether the file extention has been put in, if
% not, add!

% George O'Neill, 2024


[~,~,ext] = fileparts(file);
if isempty(ext)
    file = [file '.msh'];
end

gmsh = fopen(file, 'wt');
if (gmsh < 0)
    error('File could not be written to')
end

fprintf ( gmsh, '$MeshFormat\n' );
fprintf ( gmsh, '2.0 0 8\n' );
fprintf ( gmsh, '$EndMeshFormat\n' );

nN = size(mesh.pos,1);
fprintf ( gmsh, '$Nodes\n' );
fprintf ( gmsh, '%d\n', nN );

for ii = 1 : nN
    fprintf (gmsh, '%d', ii );
    fprintf(gmsh, ' %g %g %g\n', mesh.pos(ii,1), mesh.pos(ii,2), mesh.pos(ii,3));
end
fprintf ( gmsh, '$EndNodes\n' );

% check if this tetrahedral or hexahedral, respond accordingly;
if isfield(mesh,'tet')
    elementID = 4;
    nE = size(mesh.tet,1);
elseif isfield(mesh,'hex')
    elementID = 5;
    nE = size(mesh.hex,1);
else
    error('meshes must be tetra- or hexahedrons only!')
end

% check if field for tissue label, if not, create one, where every element
% has the tissue label 0.

if  isfield(mesh,'tissue')
    tissue = mesh.tissue-1; % Needs to be 0-indexed!
else
    tissue = zeros(nE,1);
end

fprintf (gmsh, '$Elements\n');
fprintf (gmsh, '%d\n', nE);

for ii = 1:nE
    switch elementID
        case 4
            fprintf(gmsh,'%d %d %d %d %d %d %d %d \n',...
                ii,elementID,1,tissue(ii),mesh.tet(ii,1),...
                mesh.tet(ii,2),mesh.tet(ii,3),mesh.tet(ii,4));
        case 5  
            fprintf(gmsh,'%d %d %d %d %d %d %d %d %d %d %d %d\n',...
                ii,elementID,1,tissue(ii),mesh.hex(ii,1),...
                mesh.hex(ii,2),mesh.hex(ii,3),mesh.hex(ii,4),...
                mesh.hex(ii,5),mesh.hex(ii,7),mesh.hex(ii,7),mesh.hex(ii,8));
    end
end

fprintf ( gmsh, '$EndElements\n' );
fclose ( gmsh );