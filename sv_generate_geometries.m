clearvars
close all
clc

proj_init;

triaxial = 1; % Turn on or off triaxial sensors

files = [];
files.root = proj_path;

%% Import subject from torso_tools

files.subject_mesh  = fullfile(files.root,'torso_tools','example',...
    'seated_body_registered.stl');

subject             = ft_read_headshape(files.subject_mesh);

% Downsample
p                   = [];
p.vertices          = subject.pos;
p.faces             = subject.tri;
subject             = reducepatch(p,0.6);
subject.unit        = 'mm';

%% Register the canonical throax model to subject mesh

% get the fiducuals of the torso
% - left shoulder   (pt ?)
% - right shoulder  (pt ?)
% - low spine       (pt ?)
sub_fids = [1072 -614 161
    618 -569 145
    866 -473 -330];

% Use the orginal thorax mesh to determine tranform between canonical and
% subject
S                   = [];
S.subject           = subject; % must be in units of m
S.fiducials         = sub_fids;
S.plot              = 0;
T                   = tt_register_thorax(S);

% Work out the transform to make the mesh aigned along shoulders and update
R1                  = tt_align_shoulders(T);
subject             = spm_mesh_transform(subject, R1);
T                   = R1 * T;

%% Generate sensor positions

S                   = [];
S.T                 = T;
S.subject           = subject;
S.resolution        = 30;    %% sensor spacing (in mm)
S.triaxial          = 0; % Generate a radial array
grad_1ax            = tt_generate_sensor_array(S);

S.triaxial          = 1; % Generate a triaxial array
grad_3ax            = tt_generate_sensor_array(S);

%% Try and plot the model

S                   = [];
S.subject           = subject;
S.T                 = T;
S.sensors           = grad_1ax;

tt_check_registration(S);

%% Generate the source space

S                   = [];
S.subject           = subject;
S.T                 = T;
S.width             = 1;
S.depth             = 50;
S.resolution        = 10;
S.mask              = 1;

sources = tt_generate_spine_grid(S);
sources.pos(end,:) = [];
sources.inside(end,:) = [];

scatter3(sources.pos(:,1),sources.pos(:,2),sources.pos(:,3),'y.');

%% Build the Spinal Cord and Bone Meshes based on sources

radius = 8; %% radius of enclosing cylinder in mm

% Wrap the spinal cord
Mw = tt_wrap_cord_v2(sources.pos,radius,1);
% Code to mimic Matti's alteration - first subdivide triangles
Mw = spm_mesh_refine(Mw);
% Use iso2mesh to resample (50 %) and fix
[Mw.vertices, Mw.faces] = meshresample(Mw.vertices, Mw.faces,0.5);
[Mw.vertices, Mw.faces] = meshcheckrepair(Mw.vertices, Mw.faces,'meshfix');

% Wrap the spinal cord
Mb = tt_wrap_cord_v2(sources.pos,2*radius,0);
% Code to mimic Matti's alteration - first subdivide triangles
Mb = spm_mesh_refine(Mb);
% Use iso2mesh to resample (75 %) and fix
[Mb.vertices, Mb.faces] = meshresample(Mb.vertices, Mb.faces,0.75);
[Mb.vertices, Mb.faces] = meshcheckrepair(Mb.vertices, Mb.faces,'meshfix');

%% Get the canonical meshes and transform into subject's space

[bmeshes, names] = tt_load_meshes(T);

%% Collate everything

% import and convert to meters;

mesh_wm.p = Mw.vertices / 1000;
mesh_wm.e = Mw.faces;

mesh_bone.p = Mb.vertices / 1000;
mesh_bone.e = Mb.faces;

mesh_blood.p = bmeshes{1}.vertices / 1000;
mesh_blood.e = bmeshes{1}.faces;

mesh_lungs.p = bmeshes{2}.vertices / 1000;
mesh_lungs.e = bmeshes{2}.faces;

mesh_torso.p = bmeshes{3}.vertices / 1000;
mesh_torso.e = bmeshes{3}.faces;

coils_1axis.p = grad_1ax.coilpos / 1000;
coils_1axis.n = grad_1ax.coilori;

coils_3axis.p = grad_3ax.coilpos / 1000;
coils_3axis.n = grad_3ax.coilori;

sources_cent.p = sources.pos / 1000;

sources_disk.p = disc_src / 1000; 

save(fullfile(files.root,'geometries','all_geometries.mat'),...
    'mesh_*','coils_*','sources_*');