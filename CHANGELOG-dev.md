
<a name="v0.6.0a1"></a>
## [v0.6.0a1](https://github.com/LightForm-group/cipher-parse/compare/v0.6.0a0...v0.6.0a1) - 2024.03.28

### 🐛 Bug Fixes

* use pyvista to resample VTU files onto a grid


<a name="v0.6.0a0"></a>
## [v0.6.0a0](https://github.com/LightForm-group/cipher-parse/compare/v0.5.1a4...v0.6.0a0) - 2024.03.28

### ✨ Features

* add to/from_zarr (without CIPHERGeometry support for now


<a name="v0.5.1a4"></a>
## [v0.5.1a4](https://github.com/LightForm-group/cipher-parse/compare/v0.5.1a3...v0.5.1a4) - 2024.01.08

### 🐛 Bug Fixes

* more reasonable time in notebook examples


<a name="v0.5.1a3"></a>
## [v0.5.1a3](https://github.com/LightForm-group/cipher-parse/compare/v0.5.1a2...v0.5.1a3) - 2024.01.08

### 🐛 Bug Fixes

* clear notebook outputs
* interface width is now expected within the interface definitions


<a name="v0.5.1a2"></a>
## [v0.5.1a2](https://github.com/LightForm-group/cipher-parse/compare/v0.5.1a1...v0.5.1a2) - 2023.09.05

### ✨ Features

* support separate text files for mappings

### 🐛 Bug Fixes

* DefDAP uses P=+1 for quat ops where we were assuming P=-1
* tests
* try indenting
* func args
* func args


<a name="v0.5.1a1"></a>
## [v0.5.1a1](https://github.com/LightForm-group/cipher-parse/compare/v0.5.1a0...v0.5.1a1) - 2023.07.21

### ✨ Features

* use DefDAP to compute CIPHERGeometry.misorientation_matrix - much faster


<a name="v0.5.1a0"></a>
## [v0.5.1a0](https://github.com/LightForm-group/cipher-parse/compare/v0.5.0...v0.5.1a0) - 2023.06.01

### 🐛 Bug Fixes

* cast pyvista_ndarray to numpy ndarray
* get_time_linear_subset_indices should return primitive types
* get_incremental_data if original_viz dir exists


<a name="v0.5.0"></a>
## [v0.5.0](https://github.com/LightForm-group/cipher-parse/compare/v0.5.0a3...v0.5.0) - 2023.05.15


<a name="v0.5.0a3"></a>
## [v0.5.0a3](https://github.com/LightForm-group/cipher-parse/compare/v0.5.0a2...v0.5.0a3) - 2023.05.15

### ♻ Code Refactoring

* show_interface_energies_by_misorientation to plot mobility and/or energy
* return ori range and idx from sample_from_orientations_gradient

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

### 🐛 Bug Fixes

* demo_cipher_inputs notebook
* CIPHERGeometry.from_JSON.from_JSON array casting
* index problem in CIPHERGeometry.phase_orientation setter
* CIPHERGeometry track is_periodic
* VoxelMap.get_interface_idx and VoxelMap.get_neighbour_region when not periodic

### 👷 Build changes

* update deps
* fix deps problems
* update poetry lock


<a name="v0.5.0a2"></a>
## [v0.5.0a2](https://github.com/LightForm-group/cipher-parse/compare/v0.5.0a1...v0.5.0a2) - 2023.03.03

### ✨ Features

* add ability to sample orientations from a linear orientation gradient

### 🐛 Bug Fixes

* CIPHEROutput.to_JSON and CIPHERGeometry._get_phase_orientation


<a name="v0.5.0a1"></a>
## [v0.5.0a1](https://github.com/LightForm-group/cipher-parse/compare/v0.5.0a0...v0.5.0a1) - 2023.02.15

### ♻ Code Refactoring

* separate move classes into own modules
* move CIPHERGeometry to its own module

### ✨ Features

* add CIPHEROutput.show_slice_evolution and supporting code
* identify grain boundaries in CIPHERGeometry
* add VoxelMap.get_coordinates
* CIPHEROutput.get_geometry to load an evolved geometry
* add CIPHERGeometry.get_slice and .show_slice
* allow loading CIPHERGeometry where some phases are missing
* add quiet bool arg to suppress some prints
* initial add of bin_interfaces_by_misorientation_angle

### 🐛 Bug Fixes

* _geometries assignment
* add CIPHEROutput.geometries prop
* CIPHERGeometry to/from JSON with GBs
* exception message in CIPHERGeometry init
* generate_VTI_files_from_VTU_files for 2D geometries
* pandas deprecation warning: append -> concatentate

### 👷 Build changes

* merge from remote


<a name="v0.5.0a0"></a>
## [v0.5.0a0](https://github.com/LightForm-group/cipher-parse/compare/v0.4.0...v0.5.0a0) - 2023.01.22

### ✨ Features

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

* constrain ipywidgets
* deps update
* deps
* merge branch 'aplowman/develop' of https://github.com/LightForm-group/cipher-parse into aplowman/develop
* update deps
* add notebook extras
* add optional matflow extra dep


<a name="v0.4.0"></a>
## [v0.4.0](https://github.com/LightForm-group/cipher-parse/compare/v0.3.1a3...v0.4.0) - 2022.10.26


<a name="v0.3.1a3"></a>
## [v0.3.1a3](https://github.com/LightForm-group/cipher-parse/compare/v0.3.1a2...v0.3.1a3) - 2022.10.26

### ✨ Features

* add util funcs for selecting subset of VTKs
* add func parse_cipher_stdout

### 🐛 Bug Fixes

* **build:** do not constrain poetry in release

### 👷 Build changes

* update deps


<a name="v0.3.1a2"></a>
## [v0.3.1a2](https://github.com/LightForm-group/cipher-parse/compare/v0.3.1a1...v0.3.1a2) - 2022.07.27

### 🐛 Bug Fixes

* get_misorientation_matrix orientation order


<a name="v0.3.1a1"></a>
## [v0.3.1a1](https://github.com/LightForm-group/cipher-parse/compare/v0.3.1a0...v0.3.1a1) - 2022.07.14

### ✨ Features

* add keep_arrays bool to to_JSON methods


<a name="v0.3.1a0"></a>
## [v0.3.1a0](https://github.com/LightForm-group/cipher-parse/compare/v0.3.0...v0.3.1a0) - 2022.07.13

### 🐛 Bug Fixes

* flush print in get_misorientation_matrix
* IntDefn to_JSON when metadata


<a name="v0.3.0"></a>
## [v0.3.0](https://github.com/LightForm-group/cipher-parse/compare/v0.2.3a7...v0.3.0) - 2022.07.10


<a name="v0.2.3a7"></a>
## [v0.2.3a7](https://github.com/LightForm-group/cipher-parse/compare/v0.2.3a6...v0.2.3a7) - 2022.07.10

### 🐛 Bug Fixes

* release workflow on compatible py vers


<a name="v0.2.3a6"></a>
## [v0.2.3a6](https://github.com/LightForm-group/cipher-parse/compare/v0.2.3a5...v0.2.3a6) - 2022.07.10

### ✨ Features

* add to/from_JSON to CIPHERInput

### 🐛 Bug Fixes

* PhaseTypeDefinition.to_JSON w/ orientations
* InterfaceDefinition mat/phase type tuple
* add example data

### 👷 Build changes

* test on py3.8 only
* update lock file
* downgrade h5py dep for matflow


<a name="v0.2.3a5"></a>
## [v0.2.3a5](https://github.com/LightForm-group/cipher-parse/compare/v0.2.3a4...v0.2.3a5) - 2022.07.06

### 👷 Build changes

* add itkwidgets


<a name="v0.2.3a4"></a>
## [v0.2.3a4](https://github.com/LightForm-group/cipher-parse/compare/v0.2.3a3...v0.2.3a4) - 2022.07.06

### 🐛 Bug Fixes

* add apt.txt for binder


<a name="v0.2.3a3"></a>
## [v0.2.3a3](https://github.com/LightForm-group/cipher-parse/compare/v0.2.3a2...v0.2.3a3) - 2022.07.06

### 🐛 Bug Fixes

* simplify


<a name="v0.2.3a2"></a>
## [v0.2.3a2](https://github.com/LightForm-group/cipher-parse/compare/v0.2.3a1...v0.2.3a2) - 2022.07.06


<a name="v0.2.3a1"></a>
## [v0.2.3a1](https://github.com/LightForm-group/cipher-parse/compare/v0.2.3a0...v0.2.3a1) - 2022.07.06

### 🐛 Bug Fixes

* remove file [skip ci]
* again


<a name="v0.2.3a0"></a>
## [v0.2.3a0](https://github.com/LightForm-group/cipher-parse/compare/v0.2.2...v0.2.3a0) - 2022.07.06

### 🐛 Bug Fixes

* workflow again [skip ci]


<a name="v0.2.2"></a>
## [v0.2.2](https://github.com/LightForm-group/cipher-parse/compare/v0.2.2a2...v0.2.2) - 2022.07.06


<a name="v0.2.2a2"></a>
## [v0.2.2a2](https://github.com/LightForm-group/cipher-parse/compare/v0.2.2a1...v0.2.2a2) - 2022.07.06

### 🐛 Bug Fixes

* again workflow


<a name="v0.2.2a1"></a>
## [v0.2.2a1](https://github.com/LightForm-group/cipher-parse/compare/v0.2.2a0...v0.2.2a1) - 2022.07.06

### 🐛 Bug Fixes

* workflow
* update workflow [skip ci]


<a name="v0.2.2a0"></a>
## [v0.2.2a0](https://github.com/LightForm-group/cipher-parse/compare/v0.2.1...v0.2.2a0) - 2022.07.06

### 🐛 Bug Fixes

* workflows add reqs.txt


<a name="v0.2.1"></a>
## [v0.2.1](https://github.com/LightForm-group/cipher-parse/compare/v0.2.0...v0.2.1) - 2022.07.06

### 🐛 Bug Fixes

* remove hashes from reqs.txt

### 👷 Build changes

* add poetry exported req. [skip release]


<a name="v0.2.0"></a>
## [v0.2.0](https://github.com/LightForm-group/cipher-parse/compare/v0.2.0a0...v0.2.0) - 2022.07.06


<a name="v0.2.0a0"></a>
## [v0.2.0a0](https://github.com/LightForm-group/cipher-parse/compare/v0.1.1...v0.2.0a0) - 2022.07.06

### ✨ Features

* update README


<a name="v0.1.1"></a>
## [v0.1.1](https://github.com/LightForm-group/cipher-parse/compare/v0.1.1a0...v0.1.1) - 2022.07.06


<a name="v0.1.1a0"></a>
## [v0.1.1a0](https://github.com/LightForm-group/cipher-parse/compare/v0.1.0...v0.1.1a0) - 2022.07.06

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

