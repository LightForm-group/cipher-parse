import copy
import json
from pathlib import Path
from dataclasses import dataclass
from textwrap import indent
from typing import Optional, List, Union, Tuple, Dict

import numpy as np
import h5py
from parse import parse
from ruamel.yaml import YAML
from ruamel.yaml.scalarstring import LiteralScalarString

from cipher_parse.geometry import CIPHERGeometry
from cipher_parse.interface import InterfaceDefinition
from cipher_parse.material import MaterialDefinition
from cipher_parse.utilities import set_by_path, read_shockley, grain_boundary_mobility


def compress_1D_array(arr):
    vals = []
    nums = []
    for idx, i in enumerate(arr):
        if idx == 0:
            vals.append(i)
            nums.append(1)
            continue

        if i == vals[-1]:
            nums[-1] += 1
        else:
            vals.append(i)
            nums.append(1)

    assert sum(nums) == arr.size

    return nums, vals


def compress_1D_array_string(arr, item_delim="\n"):
    out = []
    for n, v in zip(*compress_1D_array(arr)):
        out.append(f"{n} of {v}" if n > 1 else f"{v}")

    return item_delim.join(out)


def decompress_1D_array_string(arr_str, item_delim="\n"):
    out = []
    for i in arr_str.split(item_delim):
        if not i:
            continue
        if "of" in i:
            n, i = i.split("of")
            i = [int(i.strip()) for _ in range(int(n.strip()))]
        else:
            i = [int(i.strip())]
        out.extend(i)
    return np.array(out)


@dataclass
class CIPHERInput:
    geometry: CIPHERGeometry
    components: List
    outputs: List
    solution_parameters: Dict
    quiet: Optional[bool] = False

    def __post_init__(self):
        self._validate()

    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            return False
        if (
            self.components == other.components
            and self.solution_parameters == other.solution_parameters
            and self.outputs == other.outputs
            and self.geometry == other.geometry
        ):
            return True
        return False

    def _validate(self):
        check_grid_size = (
            np.array(self.solution_parameters["initblocksize"])
            * 2 ** self.solution_parameters["initrefine"]
        )
        if not np.all(check_grid_size == np.array(self.geometry.grid_size)):
            raise ValueError(
                f"`grid_size` (specifed: {self.geometry.grid_size}) must be equal to: "
                f"`initblocksize` (specified: {self.solution_parameters['initblocksize']}) "
                f"multiplied by 2 raised to the power of `initrefine` (specified: "
                f"{self.solution_parameters['initrefine']}), calculated to be: "
                f"{check_grid_size}."
            )

    def to_JSON_file(self, path):
        data = self.to_JSON()
        path = Path(path)
        with Path(path).open("wt") as fp:
            json.dump(data, fp)
        return path

    @classmethod
    def from_JSON_file(cls, path):
        with Path(path).open("rt") as fp:
            data = json.load(fp)
        return cls.from_JSON(data)

    def to_JSON(self, keep_arrays=False):
        data = {
            "geometry": self.geometry.to_JSON(keep_arrays),
            "components": self.components,
            "outputs": self.outputs,
            "solution_parameters": self.solution_parameters,
        }
        return data

    @classmethod
    def from_JSON(cls, data, quiet=True):
        data = {
            "geometry": CIPHERGeometry.from_JSON(data["geometry"], quiet=quiet),
            "components": data["components"],
            "outputs": data["outputs"],
            "solution_parameters": data["solution_parameters"],
        }
        return cls(**data, quiet=quiet)

    @classmethod
    def get_input_maps_from_files(
        cls,
        inp_file_str,
        directory,
        get_voxel_phase=True,
        get_phase_material=True,
        get_interface=True,
    ):
        directory = Path(directory)
        voxel_phase_file = Path(directory / "voxel_phase_mapping.txt")
        phase_material_file = Path(directory / "phase_material_mapping.txt")
        interface_file = Path(directory / "interface_mapping.txt")

        inp_dat = None
        voxel_phase = None
        phase_material = None
        interface_map = None
        if any(
            i.exists() for i in (voxel_phase_file, phase_material_file, interface_file)
        ):
            inp_dat = cls.read_input_YAML_string(
                file_str=inp_file_str,
                parse_interface_map=False,
            )
        if get_voxel_phase and voxel_phase_file.exists():
            with voxel_phase_file.open("rt") as fp:
                voxel_phase_str = "".join(fp.readlines())
                voxel_phase = decompress_1D_array_string(voxel_phase_str)
                voxel_phase = voxel_phase.reshape(inp_dat["grid_size"], order="F") - 1

        if get_phase_material and phase_material_file.exists():
            with phase_material_file.open("rt") as fp:
                phase_material_str = "".join(fp.readlines())
                phase_material = decompress_1D_array_string(phase_material_str) - 1

        if get_interface and interface_file.exists():
            num_phases = inp_dat["num_phases"]
            with interface_file.open("rt") as fp:
                interface_str = "".join(fp.readlines())
                interface_map = decompress_1D_array_string(interface_str)
                interface_map = interface_map.reshape((num_phases, num_phases)) - 1
                interface_map[np.tril_indices(num_phases)] = -1  # only need one half

        return (voxel_phase, phase_material, interface_map)

    @classmethod
    def from_input_YAML_file(cls, path):
        """Generate a CIPHERInput object from a CIPHER input YAML file."""

        with Path(path).open("rt") as fp:
            file_str = "".join(fp.readlines())

        (
            voxel_phase,
            phase_material,
            interface_map,
        ) = CIPHERInput.get_input_maps_from_files(
            inp_file_str=file_str,
            directory=path.parent,
        )

        return cls.from_input_YAML_str(
            file_str=file_str,
            input_map_voxel_phase=voxel_phase,
            input_map_phase_material=phase_material,
            input_map_interface=interface_map,
        )

    @classmethod
    def read_input_YAML_file(cls, path):
        with Path(path).open("rt") as fp:
            file_str = "".join(fp.readlines())

        return cls.read_input_YAML_string(file_str=file_str)

    @staticmethod
    def read_input_YAML_string(file_str, parse_interface_map=True):
        yaml = YAML(typ="safe")
        data = yaml.load(file_str)

        header = data["header"]
        grid_size = header["grid"]
        size = header["size"]
        num_phases = header["n_phases"]

        voxel_phase = None
        unique_phase_IDs = None
        if data["mappings"]["voxel_phase_mapping"] != "voxel_phase_mapping.txt":
            voxel_phase = decompress_1D_array_string(
                data["mappings"]["voxel_phase_mapping"]
            )
            voxel_phase = voxel_phase.reshape(grid_size, order="F") - 1
            unique_phase_IDs = np.unique(voxel_phase)
            assert len(unique_phase_IDs) == num_phases

        interface_map = None
        if data["mappings"]["interface_mapping"] != "interface_mapping.txt":
            if parse_interface_map:
                interface_map = decompress_1D_array_string(
                    data["mappings"]["interface_mapping"]
                )
                interface_map = interface_map.reshape((num_phases, num_phases)) - 1
                interface_map[np.tril_indices(num_phases)] = -1  # only need one half

        phase_material = None
        if data["mappings"]["phase_material_mapping"] != "phase_material_mapping.txt":
            phase_material = (
                decompress_1D_array_string(data["mappings"]["phase_material_mapping"]) - 1
            )

        return {
            "header": header,
            "grid_size": grid_size,
            "size": size,
            "num_phases": num_phases,
            "voxel_phase": voxel_phase,
            "unique_phase_IDs": unique_phase_IDs,
            "material": data["material"],
            "interface": data["interface"],
            "interface_map": interface_map,
            "phase_material": phase_material,
            "solution_parameters": data["solution_parameters"],
        }

    @classmethod
    def from_input_YAML_str(
        cls,
        file_str,
        input_map_voxel_phase=None,
        input_map_phase_material=None,
        input_map_interface=None,
        quiet=False,
    ):
        """Generate a CIPHERInput object from a CIPHER input YAML file string."""

        yaml_dat = cls.read_input_YAML_string(file_str)
        if input_map_voxel_phase is not None:
            yaml_dat["voxel_phase"] = input_map_voxel_phase
            yaml_dat["unique_phase_IDs"] = np.unique(input_map_voxel_phase)

        if input_map_phase_material is not None:
            yaml_dat["phase_material"] = input_map_phase_material

        if input_map_interface is not None:
            yaml_dat["interface_map"] = input_map_interface

        materials = [
            MaterialDefinition(
                name=name,
                properties=dict(props),
                phases=np.where(yaml_dat["phase_material"] == idx)[0],
            )
            for idx, (name, props) in enumerate(yaml_dat["material"].items())
        ]
        interfaces = []
        for idx, (int_name, props) in enumerate(yaml_dat["interface"].items()):
            phase_pairs = np.vstack(np.where(yaml_dat["interface_map"] == idx)).T
            if phase_pairs.size:
                mat_1 = materials[yaml_dat["phase_material"][phase_pairs[0, 0]]].name
                mat_2 = materials[yaml_dat["phase_material"][phase_pairs[0, 1]]].name
                type_label_part = parse(f"{mat_1}-{mat_2}{{}}", int_name)
                type_label = None
                if type_label_part:
                    type_label = type_label_part[0].lstrip("-")
                interfaces.append(
                    InterfaceDefinition(
                        properties=dict(props),
                        phase_pairs=phase_pairs,
                        materials=(mat_1, mat_2),
                        type_label=type_label,
                    )
                )

        geom = CIPHERGeometry(
            materials=materials,
            interfaces=interfaces,
            voxel_phase=yaml_dat["voxel_phase"],
            size=yaml_dat["size"],
            quiet=quiet,
        )

        attrs = {
            "geometry": geom,
            "components": yaml_dat["header"]["components"],
            "outputs": yaml_dat["header"]["outputs"],
            "solution_parameters": dict(yaml_dat["solution_parameters"]),
        }

        return cls(**attrs, quiet=quiet)

    @classmethod
    def from_voronoi(
        cls,
        grid_size,
        size,
        materials,
        interfaces,
        components,
        outputs,
        solution_parameters,
        seeds=None,
        num_phases=None,
        random_seed=None,
        is_periodic=False,
    ):
        geometry = CIPHERGeometry.from_voronoi(
            num_phases=num_phases,
            seeds=seeds,
            interfaces=interfaces,
            materials=materials,
            grid_size=grid_size,
            size=size,
            random_seed=random_seed,
            is_periodic=is_periodic,
        )

        inp = cls(
            geometry=geometry,
            components=components,
            outputs=outputs,
            solution_parameters=solution_parameters,
        )
        return inp

    @classmethod
    def from_seed_voronoi(
        cls,
        seeds,
        grid_size,
        size,
        materials,
        interfaces,
        components,
        outputs,
        solution_parameters,
        random_seed=None,
        is_periodic=False,
    ):
        return cls.from_voronoi(
            seeds=seeds,
            grid_size=grid_size,
            size=size,
            materials=materials,
            interfaces=interfaces,
            components=components,
            outputs=outputs,
            solution_parameters=solution_parameters,
            random_seed=random_seed,
            is_periodic=is_periodic,
        )

    @classmethod
    def from_random_voronoi(
        cls,
        num_phases,
        grid_size,
        size,
        materials,
        interfaces,
        components,
        outputs,
        solution_parameters,
        random_seed=None,
        is_periodic=False,
    ):
        return cls.from_voronoi(
            num_phases=num_phases,
            grid_size=grid_size,
            size=size,
            materials=materials,
            interfaces=interfaces,
            components=components,
            outputs=outputs,
            solution_parameters=solution_parameters,
            random_seed=random_seed,
            is_periodic=is_periodic,
        )

    @classmethod
    def from_voxel_phase_map(
        cls,
        voxel_phase,
        size,
        materials,
        interfaces,
        components,
        outputs,
        solution_parameters,
        random_seed=None,
    ):
        geometry = CIPHERGeometry(
            voxel_phase=voxel_phase,
            materials=materials,
            interfaces=interfaces,
            size=size,
            random_seed=random_seed,
        )
        inp = cls(
            geometry=geometry,
            components=components,
            outputs=outputs,
            solution_parameters=solution_parameters,
        )
        return inp

    @classmethod
    def from_dream3D(
        cls,
        path,
        materials,
        interfaces,
        components,
        outputs,
        solution_parameters,
        container_labels=None,
        phase_type_map=None,
    ):
        default_container_labels = {
            "SyntheticVolumeDataContainer": "SyntheticVolumeDataContainer",
            "CellData": "CellData",
            "CellEnsembleData": "CellEnsembleData",
            "FeatureIds": "FeatureIds",
            "Grain Data": "Grain Data",
            "Phases": "Phases",
            "NumFeatures": "NumFeatures",
            "BoundaryCells": "BoundaryCells",
            "NumNeighbors": "NumNeighbors",
            "NeighborList": "NeighborList",
            "SharedSurfaceAreaList": "SharedSurfaceAreaList",
            "SurfaceFeatures": "SurfaceFeatures",
            "AvgQuats": "AvgQuats",
        }
        container_labels = container_labels or {}
        container_labels = {**default_container_labels, **container_labels}

        with h5py.File(path, "r") as fp:
            voxel_phase_path = "/".join(
                (
                    "DataContainers",
                    container_labels["SyntheticVolumeDataContainer"],
                    container_labels["CellData"],
                    container_labels["FeatureIds"],
                )
            )
            phase_material_path = "/".join(
                (
                    "DataContainers",
                    container_labels["SyntheticVolumeDataContainer"],
                    container_labels["Grain Data"],
                    container_labels["Phases"],
                )
            )
            spacing_path = "/".join(
                (
                    "DataContainers",
                    container_labels["SyntheticVolumeDataContainer"],
                    "_SIMPL_GEOMETRY",
                    "SPACING",
                )
            )
            dims_path = "/".join(
                (
                    "DataContainers",
                    container_labels["SyntheticVolumeDataContainer"],
                    "_SIMPL_GEOMETRY",
                    "DIMENSIONS",
                )
            )
            material_names_path = "/".join(
                (
                    "DataContainers",
                    container_labels["SyntheticVolumeDataContainer"],
                    container_labels["CellEnsembleData"],
                    "PhaseName",
                )
            )
            grain_quats_path = "/".join(
                (
                    "DataContainers",
                    container_labels["SyntheticVolumeDataContainer"],
                    container_labels["Grain Data"],
                    container_labels["AvgQuats"],
                )
            )

            voxel_phase = fp[voxel_phase_path][()][:, :, :, 0]
            phase_material = fp[phase_material_path][()].flatten()
            voxel_phase = np.transpose(voxel_phase, axes=[2, 1, 0])
            spacing = fp[spacing_path][()]  # same as "resolution" in GUI
            dimensions = fp[dims_path][()]
            size = np.array([i * j for i, j in zip(spacing, dimensions)])
            mat_names = [i.decode("utf-8") for i in fp[material_names_path]]
            grain_quats = fp[grain_quats_path][()]

        # ignore unknown phase:
        phase_material = phase_material[1:] - 1
        grain_quats = grain_quats[1:]
        voxel_phase = voxel_phase - 1
        mat_names = mat_names[1:]

        for mat_idx, mat_name_i in enumerate(mat_names):
            phases_set = False
            if phase_type_map:
                phase_type_name = phase_type_map[mat_name_i]
            else:
                phase_type_name = mat_name_i
            for mat in materials:
                for phase_type_i in mat.phase_types:
                    if phase_type_i.name == phase_type_name:
                        phase_i_idx = np.where(phase_material == mat_idx)[0]
                        phase_type_i.phases = phase_i_idx
                        phase_type_i.orientations = grain_quats[phase_i_idx]
                        phases_set = True
                        break
                if phases_set:
                    break

            if not phases_set:
                raise ValueError(
                    f"No defined material/phase-type for Dream3D phase {mat_name_i!r}"
                )

        return cls.from_voxel_phase_map(
            voxel_phase=voxel_phase,
            size=size,
            materials=materials,
            interfaces=interfaces,
            components=components,
            outputs=outputs,
            solution_parameters=solution_parameters,
        )

    @property
    def materials(self):
        return self.geometry.materials

    @property
    def material_properties(self):
        return self.geometry.material_properties

    @property
    def interfaces(self):
        return self.geometry.interfaces

    @property
    def interface_names(self):
        return self.geometry.interface_names

    def get_header(self):
        out = {
            "grid": self.geometry.grid_size.tolist(),
            "size": self.geometry.size.tolist(),
            "n_phases": self.geometry.num_phases,
            "materials": self.geometry.material_names,
            "interfaces": self.geometry.interface_names,
            "components": self.components,
            "outputs": self.outputs,
        }
        return out

    def get_interfaces(self):
        return {i.name: i.properties for i in self.geometry.interfaces}

    def write_yaml(self, path, separate_mappings=False):
        """Write the CIPHER input YAML file.

        Parameters
        ----------
        separate_mappings
            If True, write separate text files for the mappings.

        """

        self.geometry._validate_interface_map()

        phase_mat_str = compress_1D_array_string(self.geometry.phase_material + 1) + "\n"
        vox_phase_str = (
            compress_1D_array_string(self.geometry.voxel_phase.flatten(order="F") + 1)
            + "\n"
        )
        int_str = (
            compress_1D_array_string(self.geometry.interface_map_int.flatten() + 1) + "\n"
        )

        if separate_mappings:
            phase_mat_map = "phase_material_mapping.txt"
            vox_phase_map = "voxel_phase_mapping.txt"
            int_map = "interface_mapping.txt"
            with Path(phase_mat_map).open("wt") as fh:
                fh.write(indent(phase_mat_str, "    "))

            with Path(vox_phase_map).open("wt") as fh:
                fh.write(indent(vox_phase_str, "    "))

            with Path(int_map).open("wt") as fh:
                fh.write(indent(int_str, "    "))

        else:
            phase_mat_map = LiteralScalarString(phase_mat_str)
            vox_phase_map = LiteralScalarString(vox_phase_str)
            int_map = LiteralScalarString(int_str)

        cipher_input_data = {
            "header": self.get_header(),
            "solution_parameters": dict(sorted(self.solution_parameters.items())),
            "material": {
                k: copy.deepcopy(v) for k, v in self.material_properties.items()
            },
            "interface": {k: copy.deepcopy(v) for k, v in self.get_interfaces().items()},
            "mappings": {
                "phase_material_mapping": phase_mat_map,
                "voxel_phase_mapping": vox_phase_map,
                "interface_mapping": int_map,
            },
        }

        yaml = YAML()
        path = Path(path)
        with path.open("wt", newline="\n") as fp:
            yaml.dump(cipher_input_data, fp)

        return path

    def bin_interfaces_by_misorientation_angle(
        self,
        base_interface_name,
        theta_max,
        energy_range=None,
        mobility_range=None,
        n=4,
        B=5,
        bin_width=5,
        degrees=True,
        **kwargs,
    ):
        if energy_range is None and mobility_range is None:
            raise ValueError(
                "Specify at least one of `energy_range` and `mobility_range`."
            )

        if self.geometry.misorientation_matrix is None:
            misori_matrix = self.geometry.get_misorientation_matrix(**kwargs)
        else:
            misori_matrix = self.geometry.misorientation_matrix

        min_mis, max_mis = np.min(misori_matrix), np.max(misori_matrix)
        min_range = np.floor(min_mis / bin_width) * bin_width
        max_range = np.ceil(max_mis / bin_width) * bin_width + bin_width
        misori_bins = np.linspace(
            min_range,
            max_range,
            num=int((max_range - min_range) / bin_width),
            endpoint=False,
        )
        bin_idx = np.digitize(
            misori_matrix,
            misori_bins,
            right=False,
        )
        theta = (misori_bins + (bin_width / 2))[:-1]

        if not isinstance(base_interface_name, list):
            base_interface_name = [base_interface_name]

        if energy_range is not None:
            energy = (
                read_shockley(
                    theta=theta,
                    E_max=(energy_range[1] - energy_range[0]),
                    theta_max=theta_max,
                    degrees=degrees,
                )
                + energy_range[0]
            )
        if mobility_range is not None:
            mobility = (
                grain_boundary_mobility(
                    theta=theta,
                    M_max=(mobility_range[1] - mobility_range[0]),
                    theta_max=theta_max,
                    degrees=degrees,
                    n=n,
                    B=B,
                )
                + mobility_range[0]
            )

        for int_name in base_interface_name:
            base_defn, phase_pairs = self.geometry.remove_interface(int_name)

            phase_pairs_bin_idx = bin_idx[phase_pairs[0], phase_pairs[1]]

            max_phase_pairs_fmt_len = 10

            num_existing_int_defns = len(self.geometry.interfaces)
            print("Preparing new interface defintions...")
            new_int_idx = 0
            for bin_idx_i, bin_i in enumerate(misori_bins, start=1):
                phase_pairs_bin_i_idx = np.where(phase_pairs_bin_idx == bin_idx_i)[0]
                if not phase_pairs_bin_i_idx.size:
                    continue

                else:
                    phase_pairs_bin_i = phase_pairs[:, phase_pairs_bin_i_idx].T
                    phase_pairs_bin_i_fmt = ",".join(
                        f"{i[0]}-{i[1]}" for i in phase_pairs_bin_i
                    )
                    if len(phase_pairs_bin_i_fmt) > max_phase_pairs_fmt_len:
                        phase_pairs_bin_i_fmt = (
                            phase_pairs_bin_i_fmt[: max_phase_pairs_fmt_len - 3] + "..."
                        )

                    props = copy.deepcopy(base_defn.properties)

                    if energy_range is not None:
                        new_e0 = energy[bin_idx_i - 1].item()
                        set_by_path(root=props, path=("energy", "e0"), value=new_e0)

                    if mobility_range is not None:
                        new_m0 = mobility[bin_idx_i - 1].item()
                        set_by_path(root=props, path=("mobility", "m0"), value=new_m0)

                    print(
                        f"  Adding {phase_pairs_bin_i_idx.size!r} phase pair(s) "
                        f"({phase_pairs_bin_i_fmt}) to bin {bin_idx_i} with mid value: "
                        f"{theta[bin_idx_i - 1]!r}."
                    )

                    new_type_lab = str(new_int_idx)
                    if base_defn.type_label:
                        new_type_lab = f"{base_defn.type_label}-{new_type_lab}"

                    new_int = InterfaceDefinition(
                        phase_types=base_defn.phase_types,
                        type_label=new_type_lab,
                        properties=props,
                        phase_pairs=phase_pairs_bin_i.tolist(),
                    )
                    self.geometry.interfaces.append(new_int)
                    self.geometry._modify_interface_map(
                        phase_A=phase_pairs_bin_i[:, 0],
                        phase_B=phase_pairs_bin_i[:, 1],
                        interface_idx=(num_existing_int_defns + new_int_idx),
                    )
                    new_int_idx += 1

        print("done!")
        self.geometry._check_interface_phase_pairs()
        self.geometry._validate_interfaces()
        self.geometry._validate_interface_map()

    def apply_interface_property(
        self,
        base_interface_name,
        property_name,
        property_values,
        additional_metadata=None,
        bin_edges=None,
    ):
        """Expand a base interface into multiple interfaces, by assigning the specified
        property (e.g. GB energy) from a symmetric matrix of property values

        Parameters
        ----------
        base_interface_name : str
        property_name : tuple of str
        property_values : ndarray of shape (N_phases, N_phases)
            N_phases it the total number of phases in the geometry.
        bin_edges : ndarray of float, optional
            If specified, bin property values such that multiple phase-pairs are
            represented by the same interface definition. This uses `np.digitize`. The
            values used for each bin will be the mid-points between bin edges, where a
            given mid-point is larger than its associated edge.

        """

        if not isinstance(property_name, list):
            property_name = [property_name]

        if not isinstance(property_values, list):
            property_values = [property_values]

        if not isinstance(bin_edges, list):
            bin_edges = [bin_edges]

        base_defn, phase_pairs = self.geometry.remove_interface(base_interface_name)
        new_vals_all = property_values[0][phase_pairs[0], phase_pairs[1]]

        new_interfaces_data = []
        if bin_edges[0] is not None:
            bin_idx = np.digitize(new_vals_all, bin_edges[0])
            all_pp_idx_i = []
            for idx, bin_i in enumerate(bin_edges[0]):
                pp_idx_i = np.where(bin_idx == idx + 1)[0]
                all_pp_idx_i.extend(pp_idx_i.tolist())
                if pp_idx_i.size:
                    if idx < len(bin_edges[0]) - 1:
                        value = (bin_i + bin_edges[0][idx + 1]) / 2
                    else:
                        value = bin_i
                    print(
                        f"Adding {pp_idx_i.size!r} phase pair(s) to {property_name!r} bin "
                        f"{idx + 1} with edge value: {bin_i!r} and centre: {value!r}."
                    )
                    new_interfaces_data.append(
                        {
                            "phase_pairs": phase_pairs.T[pp_idx_i],
                            "values": [value],
                            "bin_idx": idx,
                        }
                    )

            for name, vals, edges in zip(
                property_name[1:], property_values[1:], bin_edges[1:]
            ):
                for idx, new_int_dat in enumerate(new_interfaces_data):
                    bin_idx = new_int_dat["bin_idx"]
                    bin_i = edges[bin_idx]
                    if bin_idx < len(edges) - 1:
                        value = (bin_i + edges[bin_idx + 1]) / 2
                    else:
                        value = bin_i
                    new_interfaces_data[idx]["values"].append(value)

            miss_phase_pairs = set(np.arange(phase_pairs.shape[1])) - set(all_pp_idx_i)
            if miss_phase_pairs:
                miss_prop_vals = [
                    property_values[phase_pairs[0, i], phase_pairs[1, i]]
                    for i in miss_phase_pairs
                ]
                missing_dat = dict(
                    zip(
                        (tuple(phase_pairs[:, i]) for i in miss_phase_pairs),
                        miss_prop_vals,
                    )
                )
                raise RuntimeError(
                    f"Not all phase pairs have been added to a property value bin. The "
                    f"following {len(missing_dat)}/{phase_pairs.shape[1]} phase pairs (and "
                    f"property values) are missing: {missing_dat}."
                )
        else:
            print(
                f"Adding a new interface for each of {phase_pairs.shape[1]} phase pairs."
            )
            new_interfaces_data = [
                {
                    "phase_pairs": np.array([pp]),
                    "values": [i[pp_idx] for i in new_vals_all],
                }
                for pp_idx, pp in enumerate(phase_pairs.T)
            ]

        num_existing_int_defns = len(self.geometry.interfaces)
        print("Preparing new interface defintions...", end="")
        for idx, i in enumerate(new_interfaces_data):
            props = copy.deepcopy(base_defn.properties)
            for name, val in zip(property_name, i["values"]):
                new_value = val.item()  #  convert from numpy to native
                set_by_path(root=props, path=name, value=new_value)

            new_type_lab = str(idx)
            if base_defn.type_label:
                new_type_lab = f"{base_defn.type_label}-{new_type_lab}"

            metadata = {}
            if additional_metadata:
                for k, v in additional_metadata.items():
                    metadata[k] = v[i["phase_pairs"][:, 0], i["phase_pairs"][:, 1]]

            new_int = InterfaceDefinition(
                phase_types=base_defn.phase_types,
                type_label=new_type_lab,
                properties=props,
                phase_pairs=i["phase_pairs"].tolist(),
                metadata=metadata,
            )
            self.geometry.interfaces.append(new_int)
            self.geometry._modify_interface_map(
                phase_A=i["phase_pairs"][:, 0],
                phase_B=i["phase_pairs"][:, 1],
                interface_idx=(num_existing_int_defns + idx),
            )

        print("done!")
        self.geometry._check_interface_phase_pairs()
        self.geometry._validate_interfaces()
        self.geometry._validate_interface_map()
