
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

