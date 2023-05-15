
<a name="v0.5.0"></a>
## [v0.5.0](https://github.com/LightForm-group/cipher-parse/compare/v0.4.0...v0.5.0) - 2023.05.15

### ♻ Code Refactoring

* show_interface_energies_by_misorientation to plot mobility and/or energy
* return ori range and idx from sample_from_orientations_gradient
* separate move classes into own modules
* move CIPHERGeometry to its own module

### ✨ Features

* add IPF_z data label type to CIPHERGeometry.get_slice
* add discrete_colours arg to CIPHERGeometry.show_slice
* add interface data label type to CIPHERGeometry.get_slice
* add CIPHERGeometry.show_relative_misorientation
* plot energy and/or mobility in generate_interface_energies_plot
* add VoxelMap.get_interface_voxels
* add CIPHEROutput.get_average_radius_evolution
* add discrete_colours arg to CIPHEROutput.show_slice_evolution
* support multiple interface_binning base interfaces
* allow energy or/and mobility ranges in bin_interfaces_by_misorientation_angle
* add CIPHERGeometry get/show_interface_energies_by_misorientation
* add CIPHERGeometry get/show_interface_energies_by_misorientation
* add increment and incremental_data_idx to CIPHERGeometry
* show sim time on show_phase_size_dist_evolution x-axis
* allow passing CIPHERInput to CIPHEROutput init
* add ability to sample orientations from a linear orientation gradient
* add CIPHEROutput.show_slice_evolution and supporting code
* identify grain boundaries in CIPHERGeometry
* add VoxelMap.get_coordinates
* CIPHEROutput.get_geometry to load an evolved geometry
* add CIPHERGeometry.get_slice and .show_slice
* allow loading CIPHERGeometry where some phases are missing
* add quiet bool arg to suppress some prints
* initial add of bin_interfaces_by_misorientation_angle
* add grain_boundary_mobility empirical form
* support saving outputs via time increment
* save misori matrix
* add compare_phase_size_dist_evolution
* add bin_size to show_phase_size_dist_evo
* add max_increments opt to show_phase_size_dist_evolution
* add `read_input_YAML_string`
* add option to show size dist as prob
* add CIPHERGeometry.write_VTK; fix [#35](https://github.com/LightForm-group/cipher-parse/issues/35)
* select orientations when assigning phases
* add options to CIPHEROutput
* add derived outputs
* add show_phase_size_dist_evolution
* add CIPHEROutput
* read YAML file back in
* set int phase_pairs on _get_interface_map
* add __eq__ to CIPHERInput classes
* add is_periodic to voronoi gen

### 🐛 Bug Fixes

* demo_cipher_inputs notebook
* CIPHERGeometry.from_JSON.from_JSON array casting
* index problem in CIPHERGeometry.phase_orientation setter
* CIPHERGeometry track is_periodic
* VoxelMap.get_interface_idx and VoxelMap.get_neighbour_region when not periodic
* CIPHEROutput.to_JSON and CIPHERGeometry._get_phase_orientation
* _geometries assignment
* add CIPHEROutput.geometries prop
* CIPHERGeometry to/from JSON with GBs
* exception message in CIPHERGeometry init
* generate_VTI_files_from_VTU_files for 2D geometries
* pandas deprecation warning: append -> concatentate
* CIPHERGeometry.from_JSON
* failing round trip test
* parse type label in CIPHERInput from YAML
* glob specificity for VTK search
* better output options naming
* get_subset_indices
* more robust VTK subset selector
* CIPHERInput from YAML with unused interface
* decompress array empty line bug

### 👷 Build changes

* update deps
* fix deps problems
* update poetry lock
* merge from remote
* constrain ipywidgets
* deps update
* deps
* merge branch 'aplowman/develop' of https://github.com/LightForm-group/cipher-parse into aplowman/develop
* update deps
* add notebook extras
* add optional matflow extra dep


<a name="v0.4.0"></a>
## [v0.4.0](https://github.com/LightForm-group/cipher-parse/compare/v0.3.0...v0.4.0) - 2022.10.26

### ✨ Features

* add util funcs for selecting subset of VTKs
* add func parse_cipher_stdout
* add keep_arrays bool to to_JSON methods

### 🐛 Bug Fixes

* get_misorientation_matrix orientation order
* flush print in get_misorientation_matrix
* IntDefn to_JSON when metadata
* **build:** do not constrain poetry in release

### 👷 Build changes

* update deps


<a name="v0.3.0"></a>
## [v0.3.0](https://github.com/LightForm-group/cipher-parse/compare/v0.2.2...v0.3.0) - 2022.07.10

### ✨ Features

* add to/from_JSON to CIPHERInput

### 🐛 Bug Fixes

* release workflow on compatible py vers
* PhaseTypeDefinition.to_JSON w/ orientations
* InterfaceDefinition mat/phase type tuple
* add example data
* add apt.txt for binder
* simplify
* remove file [skip ci]
* again
* workflow again [skip ci]

### 👷 Build changes

* test on py3.8 only
* update lock file
* downgrade h5py dep for matflow
* add itkwidgets


<a name="v0.2.2"></a>
## [v0.2.2](https://github.com/LightForm-group/cipher-parse/compare/v0.2.1...v0.2.2) - 2022.07.06

### 🐛 Bug Fixes

* again workflow
* workflow
* update workflow [skip ci]
* workflows add reqs.txt


<a name="v0.2.1"></a>
## [v0.2.1](https://github.com/LightForm-group/cipher-parse/compare/v0.2.0...v0.2.1) - 2022.07.06

### 🐛 Bug Fixes

* remove hashes from reqs.txt

### 👷 Build changes

* add poetry exported req. [skip release]


<a name="v0.2.0"></a>
## [v0.2.0](https://github.com/LightForm-group/cipher-parse/compare/v0.1.1...v0.2.0) - 2022.07.06

### ✨ Features

* update README


<a name="v0.1.1"></a>
## [v0.1.1](https://github.com/LightForm-group/cipher-parse/compare/v0.1.0...v0.1.1) - 2022.07.06

### 🐛 Bug Fixes

* test release
* revert vers bump
* release workflow
* test release


<a name="v0.1.0"></a>
## v0.1.0 - 2022.07.06

### ✨ Features

* add core code
* first commit

### 👷 Build changes

* change test py vers

