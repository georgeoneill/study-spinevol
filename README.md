# spinevol
 
Supporting code to produce the plots from the paper [Volume conductor models for magnetospinography](https://www.biorxiv.org/content/10.1101/2024.11.04.621905v1)

>:hammer:**Note:** This repository is under development whilst the manuscript is subject to peer review, an archived version of the code and dependencies can be found on [Zenodo](https://zenodo.org/records/14883494) once the paper has been accepted (and for any urgent fixes).

>:warning:**Warning**: There are two current issues with this repository 
**1)** The functions used to plot the spinal topoplots in the paper are currently unavailable due to licensing issues: a license friendly workaround is under construction.
**2)** Functionality to call the FEM with non-windows operating systems is not supported, the binaries exist [here](https://github.com/brainstorm-tools/bst-duneuro/tree/master/bin) but no validation that this local code works with them has been tested yet.

## Usage

#### Geometry generation
If you want to generate the source space, sensor layouts and boundaries from scratch, please run `sv_generate_geometries.m`. These will be very similar, but not completely identical to the boundaries used in the paper (due to some stochasticity in the code). To use those, please load in `geometries/all_geometries_v2.mat`.

#### Forward Model Solving
We provide two scripts to solve the forward models, `sv_make_lead_fields_central.m` for the sources which are in the centre of the spinal cord `sv_make_lead_fields_disk.m` sources in a ring away from the centre.

#### Analysis replication from paper

Once the forward models are solved, the assets from the figures in the manuscript 

Figure 2: `sv_make_topoplots.m`
Figure 3: `sv_compare_fields_central.m`
Figure 4: `sv_pca_fields.m`
Figure 5: `sv_disk_source_comps.m`

## Dependencies

The code is built upon existing open-source repositories and are included as submodule for convenience. To install including submodules use

```bash
git clone --recurse-submodules github.com/georgeoneill/study-spinevol
```

Drop the `--recurse-submodules` flag to download only the local code. 


#### Modules

If you want to hunt out the modules yourself, they are:

- [hbf_lc_p](https://github.com/MattiStenroos/hbf_lc_p): Helsinki BEM Framework LCISA solver for MEG/EEG.
- [Torso Tools](https://github.com/fil-opmeg/torso_tools): Generation of spinal source, sensor and conductive layouts, and functions to call hbf_lc_p.
- [SPM](https://github.com/spm/spm): Neuroimaging analysis software, containing mesh processing functionality and functions calls to solve forward models.
- [ISO2MESH](https://github.com/fangq/iso2mesh): Advanced meshing tools.
- [bst-duneuro](https://github.com/brainstorm-tools/bst-duneuro): compiled DuNeuro code for finite element analysis of neuroimaging models.
- [gramm](https://github.com/piermorel/gramm): A plotting library.
