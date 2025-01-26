function L = fem_calc_fwds(S)

% Helper script to run DuNeuro on a given mesh.
% Add duneuro and bst functions to path

% George O'Neill, 2024


if ~isfield(S,'dir'), S.dir = tempname; end
if ~isfield(S,'grad'), error('Please provide fieldtrip grad struct!'); end
if ~isfield(S,'src'), error('Please provide fieldtrip src struct!'); end
if ~isfield(S,'mesh'), error('Please provide fieldtrip tet/hex struct!'); end
if ~isfield(S,'cond'), error('Please provide conductivities!'); end
if ~isfield(S,'binpath'), S.binpath = []; end

% check results dir exist
if ~exist(S.dir,'dir')
    mkdir(S.dir);
end

TmpDir = S.dir;
dnModality = 'meg';

cfg = duneuro_defaults();
cfg.FemCond = S.cond;
cfg.BstSaveTransfer = 1;


% write out mesh
% Check if mesh type

if isstruct(S.mesh)
    if isfield(S.mesh,'tet')
        ElementType = 'tetrahedron';
    elseif isfield(S.mesh,'hex')
        ElementType = 'hexahedron';
    else
        error('please specify if mesh is tetrahedral or hexahedral')
    end
else
    error('mesh must be in fieldtrip structure format!')
end
MeshFile = fullfile(TmpDir,'volume_conductor.msh');
fem_write_gmsh(MeshFile,S.mesh);


% write out dipoles
DipoleFile = fullfile(TmpDir,'dipole_model.txt');
dipoles = [kron(S.src.pos, ones(3,1)), kron(ones(size(S.src.pos,1), 1), eye(3))];
fid = fopen(DipoleFile, 'wt+');
fprintf(fid, '%d %d %d %d %d %d \n', dipoles');
fclose(fid);

% write out coil positions
CoilFile = fullfile(TmpDir, 'coil_model.txt');
CoilsLoc = S.grad.coilpos;
fid = fopen(CoilFile, 'wt+');
fprintf(fid, '%d %d %d  \n', CoilsLoc');
fclose(fid);

% write out coil orientations
ProjFile = fullfile(TmpDir, 'projection_model.txt');
CoilsOrient = S.grad.coilori;
fid = fopen(ProjFile, 'wt+');
fprintf(fid, '%d %d %d  \n', CoilsOrient');
fclose(fid);

% write out conductivity data
CondFile = fullfile(TmpDir, 'conductivity_model.con');
fid = fopen(CondFile, 'w');
fprintf(fid, '%d\t', cfg.FemCond);
fclose(fid);

%% ===== WRITE MINI FILE =====
% Open the mini file
IniFile = fullfile(TmpDir, 'duneuro_minifile.mini');
fid = fopen(IniFile, 'wt+');
% General setting
fprintf(fid, '__name = %s\n\n', IniFile);
if strcmp(cfg.SolverType, 'cg')
    fprintf(fid, 'type = %s\n', cfg.FemType);
end
fprintf(fid, 'element_type = %s\n', ElementType);
fprintf(fid, 'solver_type = %s\n', cfg.SolverType);
fprintf(fid, 'geometry_adapted = %s\n', bool2str(cfg.GeometryAdapted));
fprintf(fid, 'tolerance = %d\n', cfg.Tolerance);
if strcmp(dnModality, 'eeg') || strcmp(dnModality, 'meeg')
    fprintf(fid, '[electrodes]\n');
    fprintf(fid, 'filename = %s\n', fullfile(TmpDir, ElecFile));
    fprintf(fid, 'type = %s\n', cfg.ElecType);
    fprintf(fid, 'codims = %s\n', '3');
end
% [meg]
if strcmp(dnModality, 'meg') || strcmp(dnModality, 'meeg')
    fprintf(fid, '[meg]\n');
    fprintf(fid, 'intorderadd = %d\n', cfg.MegIntorderadd);
    fprintf(fid, 'type = %s\n', cfg.MegType);
    fprintf(fid, 'cache.enable = %s\n',bool2str(cfg.EnableCacheMemory) );
    % [coils]
    fprintf(fid, '[coils]\n');
    fprintf(fid, 'filename = %s\n', CoilFile);
    % [projections]
    fprintf(fid, '[projections]\n');
    fprintf(fid, 'filename = %s\n', ProjFile);
end
% [dipoles]
fprintf(fid, '[dipoles]\n');
fprintf(fid, 'filename = %s\n', DipoleFile);
% [volume_conductor.grid]
fprintf(fid, '[volume_conductor.grid]\n');
fprintf(fid, 'filename = %s\n', MeshFile);
% [volume_conductor.tensors]
fprintf(fid, '[volume_conductor.tensors]\n');
fprintf(fid, 'filename = %s\n', CondFile);
% [solver]
fprintf(fid, '[solver]\n');
fprintf(fid, 'solver_type = %s\n', cfg.SolvSolverType);
fprintf(fid, 'preconditioner_type = %s\n', cfg.SolvPrecond);
if strcmp(cfg.SolverType, 'cg')
    fprintf(fid, 'cg_smoother_type = %s\n', cfg.SolvSmootherType);
end
fprintf(fid, 'intorderadd = %d\n', cfg.SolvIntorderadd);
% Discontinuous Galerkin
if strcmp(cfg.SolverType, 'dg')
    fprintf(fid, 'dg_smoother_type = %s\n', cfg.DgSmootherType);
    fprintf(fid, 'scheme = %s\n', cfg.DgScheme);
    fprintf(fid, 'penalty = %d\n', cfg.DgPenalty);
    fprintf(fid, 'edge_norm_type = %s\n', cfg.DgEdgeNormType);
    fprintf(fid, 'weights = %s\n', bool2str(cfg.DgWeights));
    fprintf(fid, 'reduction = %s\n', bool2str(cfg.DgReduction));
end
% [solution]
fprintf(fid, '[solution]\n');
fprintf(fid, 'post_process = %s\n', bool2str(cfg.SolPostProcess)); % true/false
fprintf(fid, 'subtract_mean = %s\n', bool2str(cfg.SolSubstractMean)); % boolean
% [solution.solver]
fprintf(fid, '[solution.solver]\n');
fprintf(fid, 'reduction = %d\n', cfg.SolSolverReduction);
% [solution.source_model]
fprintf(fid, '[solution.source_model]\n');
fprintf(fid, 'type = %s\n', cfg.SrcModel);
fprintf(fid, 'intorderadd = %d\n', cfg.SrcIntorderadd);
fprintf(fid, 'intorderadd_lb = %d\n', cfg.SrcIntorderadd_lb);
fprintf(fid, 'numberOfMoments = %d\n', cfg.SrcNbMoments);
fprintf(fid, 'referenceLength = %d\n', cfg.SrcRefLen);
fprintf(fid, 'weightingExponent = %d\n', cfg.SrcWeightExp);
fprintf(fid, 'relaxationFactor = %e\n', 10^(-cfg.SrcRelaxFactor));
fprintf(fid, 'mixedMoments = %s\n', bool2str(cfg.SrcMixedMoments));
fprintf(fid, 'restrict = %s\n', bool2str(cfg.SrcRestrict));
fprintf(fid, 'initialization = %s\n', cfg.SrcInit);
% [brainstorm]
fprintf(fid, '[brainstorm]\n');
fprintf(fid, 'modality = %s\n', dnModality);
fprintf(fid, 'output_folder = %s\n', [TmpDir, filesep]);
fprintf(fid, 'save_eeg_transfer_file = %s\n', bool2str(cfg.BstSaveTransfer));
fprintf(fid, 'save_meg_transfer_file = %s\n', bool2str(cfg.BstSaveTransfer));
fprintf(fid, 'save_meeg_transfer_file = %s\n', bool2str(cfg.BstSaveTransfer));
fprintf(fid, 'eeg_transfer_filename = %s\n', cfg.BstEegTransferFile);
fprintf(fid, 'meg_transfer_filename = %s\n', cfg.BstMegTransferFile);
fprintf(fid, 'eeg_leadfield_filename = %s\n', cfg.BstEegLfFile);
fprintf(fid, 'meg_leadfield_filename = %s\n', cfg.BstMegLfFile);
% Close file
fclose(fid);

%% Run Binary

if isempty(S.binpath)
    [p,~,~] = fileparts(mfilename('fullpath'));
    p = fullfile(p,'private');
else
    p = S.binpath;
end

if ispc
    exepath = fullfile(p,'bst_duneuro_meeg_win64.exe');
    cmd = ['"',exepath,'"',' "',IniFile,'"'];
else
    error('Mac OS and Linux support coming soon');
end

status = system(cmd);
if status ~= 0
    error('unknown snafu with DuNeuro, see command line')
end

%% Finish off lead fields

GainMeg = in_duneuro_bin(fullfile(TmpDir, cfg.BstMegLfFile))';

% === POST-PROCESS MEG LEADFIELD ===
% Compute the total magnetic field
dipoles_pos_orie = [kron(S.src.pos,ones(3,1)), kron(ones(length(S.src.pos),1), eye(3))];

% a- Compute the MEG Primary Magnetic B-field analytically (formula of Sarvas)
dip_pos = dipoles_pos_orie(:,1:3);
dip_mom = dipoles_pos_orie(:,4:6);
Bp = zeros(size(S.grad.coilpos,1), size(dip_pos,1));
for i = 1:size(S.grad.coilpos,1)
    for j = 1 : size(dip_pos,1)
        R = S.grad.coilpos(i,:);
        R_0 = dip_pos(j,:);
        A = R - R_0;
        a = norm(A);
        aa = A./(a^3);
        BpHelp = cross(dip_mom(j,:),aa);
        Bp(i,j) = BpHelp * S.grad.coilori(i, :)'; % projection of the primary B-field along the coil orientations
    end
end

% b- The total magnetic field B = Bp + Bs;
%  full B-field
Bs = GainMeg;
mu = 4*pi*1e-7; % check the value of the units maybe it needs to be mu = 4*pi*1e-7
Bfull = (mu/(4*pi)) * (Bp - Bs);

% Correct for balancing (e.g. AMM etc).
L = S.grad.tra * Bfull;

end

function str = bool2str(bool)
if bool
    str = 'true';
else
    str = 'false';
end
end
