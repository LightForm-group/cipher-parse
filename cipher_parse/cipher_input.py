import copy
import json
from pathlib import Path
from dataclasses import dataclass
from typing import Optional, List, Union, Tuple, Dict

import numpy as np
import h5py
from parse import parse
from ruamel.yaml import YAML
from ruamel.yaml.scalarstring import LiteralScalarString

from cipher_parse.geometry import CIPHERGeometry
from cipher_parse.errors import (
    MaterialPhaseTypeFractionError,
    MaterialPhaseTypeLabelError,
    MaterialPhaseTypePhasesMissingError,
)
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


class InterfaceDefinition:
    """
    Attributes
    ----------
    materials :
        Between which named materials this interface applies.  Specify this or `phase_types`.
    phase_types :
        Between which named phase types this interface applies. Specify this or `materials`.
    type_label :
        To distinguish between multiple interfaces that all apply between the same pair of
        materials
    phase_pairs :
        List of phase pair indices that should have this interface type (for manual
        specification). Can be specified as an (N, 2) array.
    """

    def __init__(
        self,
        properties: Dict,
        materials: Optional[Union[List[str], Tuple[str]]] = None,
        phase_types: Optional[Union[List[str], Tuple[str]]] = None,
        type_label: Optional[str] = None,
        type_fraction: Optional[float] = None,
        phase_pairs: Optional[np.ndarray] = None,
        metadata: Optional[Dict] = None,
    ):
        self._is_phase_pairs_set = False
        self.index = None  # assigned by parent CIPHERGeometry

        self.properties = properties
        self.materials = tuple(materials) if materials else None
        self.phase_types = tuple(phase_types) if phase_types else None
        self.type_label = type_label
        self.type_fraction = type_fraction
        self.phase_pairs = phase_pairs
        self.metadata = metadata

        self._validate()

    def __eq__(self, other):
        # note we don't check type_fraction, should we?
        if not isinstance(other, self.__class__):
            return False
        if (
            self.type_label == other.type_label
            and sorted(self.phase_types) == sorted(other.phase_types)
            and self.properties == other.properties
            and np.all(self.phase_pairs == other.phase_pairs)
        ):
            return True
        return False

    def to_JSON(self, keep_arrays=False):
        data = {
            "properties": self.properties,
            "phase_types": list(self.phase_types),
            "type_label": self.type_label,
            "type_fraction": self.type_fraction,
            "phase_pairs": self.phase_pairs if self.is_phase_pairs_set else None,
            "metadata": {k: v for k, v in (self.metadata or {}).items()} or None,
        }
        if not keep_arrays:
            if self.is_phase_pairs_set:
                data["phase_pairs"] = data["phase_pairs"].tolist()
            if self.metadata:
                data["metadata"] = {k: v.tolist() for k, v in data["metadata"].items()}
        return data

    @classmethod
    def from_JSON(cls, data):
        data = {
            "properties": data["properties"],
            "phase_types": tuple(data["phase_types"]),
            "type_label": data["type_label"],
            "type_fraction": data["type_fraction"],
            "phase_pairs": np.array(data["phase_pairs"])
            if data["phase_pairs"] is not None
            else None,
            "metadata": (
                {k: np.array(v) for k, v in data["metadata"].items()}
                if data["metadata"]
                else None
            ),
        }
        return cls(**data)

    @property
    def is_phase_pairs_set(self):
        return self._is_phase_pairs_set

    @property
    def name(self):
        return self.get_name(self.phase_types, self.type_label)

    @property
    def phase_pairs(self):
        return self._phase_pairs

    @phase_pairs.setter
    def phase_pairs(self, phase_pairs):

        if phase_pairs is not None:
            self._is_phase_pairs_set = True

        if phase_pairs is None or len(phase_pairs) == 0:
            phase_pairs = np.array([]).reshape((0, 2))
        else:
            phase_pairs = np.asarray(phase_pairs)

        if phase_pairs.shape[1] != 2:
            raise ValueError(
                f"phase_pairs should be specified as an (N, 2) array or a list of "
                f"two-element lists, but has shape: {phase_pairs.shape}."
            )

        # sort so first index is smaller:
        phase_pairs = np.sort(phase_pairs, axis=1)

        # sort by first phase index, then by second phase-idx:
        srt = np.lexsort(phase_pairs.T[::-1])
        phase_pairs = phase_pairs[srt]

        self._phase_pairs = phase_pairs

    @property
    def num_phase_pairs(self):
        return self.phase_pairs.shape[0]

    @property
    def metadata(self):
        return self._metadata

    @metadata.setter
    def metadata(self, metadata):
        if metadata is not None:
            for k, v in metadata.items():
                if len(v) != self.num_phase_pairs:
                    raise ValueError(
                        f"Item {k!r} in the `metadata` dict must have length equal to the "
                        f"number of phase pairs ({self.num_phase_pairs}) but has length: "
                        f"{len(v)}."
                    )
        self._metadata = metadata

    @staticmethod
    def get_name(phase_types, type_label):
        return (
            f"{phase_types[0]}-{phase_types[1]}"
            f"{f'-{type_label}' if type_label else ''}"
        )

    def _validate(self):
        if self.materials:
            if self.phase_types:
                raise ValueError("Specify exactly one of `materials` and `phase_types`.")
            self.phase_types = copy.copy(self.materials)

        elif not self.phase_types:
            raise ValueError("Specify exactly one of `materials` and `phase_types`.")

        if self.type_fraction is not None and self.phase_pairs.size:
            raise ValueError("Specify either `type_fraction` or `phase_pairs`.")


class MaterialDefinition:
    """Class to represent a material within a CIPHER simulation."""

    def __init__(
        self,
        name,
        properties,
        phase_types=None,
        target_volume_fraction=None,
        phases=None,
    ):

        self.name = name
        self.properties = properties
        self.target_volume_fraction = target_volume_fraction
        self._geometry = None

        if target_volume_fraction is not None and phases is not None:
            raise ValueError(
                f"Cannot specify both `target_volume_fraction` and `phases` for material "
                f"{self.name!r}."
            )  # TODO: test raise

        if target_volume_fraction is not None:
            if target_volume_fraction == 0.0 or target_volume_fraction > 1.0:
                raise ValueError(
                    f"Target volume fraction must be greater than zero and less than or "
                    f"equal to one, but specified value for material {self.name!r} was "
                    f"{target_volume_fraction!r}."
                )  # TODO: test raise

        if phases is not None:
            for i in phase_types or []:
                if i.phases is not None:
                    raise ValueError(
                        f"Cannot specify `phases` in any of the phase type definitions if "
                        f"`phases` is also specified in the material definition."
                    )  # TODO: test raise
        else:
            if phase_types:
                is_phases_given = [i.phases is not None for i in phase_types]
                if any(is_phases_given) and sum(is_phases_given) != len(phase_types):
                    raise MaterialPhaseTypePhasesMissingError(
                        f"If specifying `phases` for a phase type for material "
                        f"{self.name!r}, `phases` must be specified for all phase types."
                    )

        if phase_types is None:
            phase_types = [PhaseTypeDefinition(phases=phases)]

        if len(phase_types) > 1:
            pt_labels = [i.type_label for i in phase_types]
            if len(set(pt_labels)) < len(pt_labels):
                raise MaterialPhaseTypeLabelError(
                    f"Phase types belonging to the same material ({self.name!r}) must have "
                    f"distinct `type_label`s."
                )

        self.phase_types = phase_types

        if self.target_volume_fraction is not None:
            if self.phases is not None:
                raise ValueError(
                    f"Cannot specify both `target_volume_fraction` and `phases` for "
                    f"material {self.name!r}."
                )  # TODO: test raise

        is_type_frac = [i.target_type_fraction is not None for i in phase_types]
        if phase_types[0].phases is None:
            num_unassigned_vol = self.num_phase_types - sum(is_type_frac)
            assigned_vol = sum(i or 0.0 for i in self.target_phase_type_fractions)
            if num_unassigned_vol:
                frac = (1.0 - assigned_vol) / num_unassigned_vol
                if frac <= 0.0:
                    raise MaterialPhaseTypeFractionError(
                        f"All phase type target volume fractions must sum to one, but "
                        f"assigned target volume fractions sum to {assigned_vol} with "
                        f"{num_unassigned_vol} outstanding unassigned phase type volume "
                        f"fraction(s)."
                    )
            for i in self.phase_types:
                if i.target_type_fraction is None:
                    i.target_type_fraction = frac

            assigned_vol = sum(self.target_phase_type_fractions)
            if not np.isclose(assigned_vol, 1.0):
                raise MaterialPhaseTypeFractionError(
                    f"All phase type target type fractions must sum to one, but target "
                    f"type fractions sum to {assigned_vol}."
                )

        for i in self.phase_types:
            i._material = self

    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            return False
        if (
            self.name == other.name
            and self.properties == other.properties
            and np.all(self.phases == other.phases)
        ):
            return True
        return False

    def to_JSON(self, keep_arrays=False):
        data = {
            "name": self.name,
            "properties": self.properties,
            "phase_types": [i.to_JSON(keep_arrays) for i in self.phase_types],
        }
        return data

    @classmethod
    def from_JSON(cls, data):
        data = {
            "name": data["name"],
            "properties": data["properties"],
            "phase_types": [
                PhaseTypeDefinition.from_JSON(i) for i in data["phase_types"]
            ],
        }
        return cls(**data)

    @property
    def geometry(self):
        return self._geometry

    @property
    def num_phase_types(self):
        return len(self.phase_types)

    @property
    def target_phase_type_fractions(self):
        return [i.target_type_fraction for i in self.phase_types]

    @property
    def phases(self):
        try:
            return np.concatenate([i.phases for i in self.phase_types])
        except ValueError:
            # phases not yet assigned
            return None

    @property
    def index(self):
        """Get the index within the geometry materials list."""
        return self.geometry.materials.index(self)

    @property
    def phase_type_fractions(self):
        """Get the actual type volume (voxel) fractions within the material."""
        phase_type_fractions = []
        for i in self.phase_types:
            num_mat_voxels = self.geometry.material_num_voxels[self.index]
            pt_num_voxels = np.sum(self.geometry.phase_num_voxels[i.phases])
            phase_type_fractions.append(pt_num_voxels / num_mat_voxels)
        return np.array(phase_type_fractions)

    def assign_phases(self, phases, random_seed=None):
        """Assign given phase indices to phase types according to target_type_fractions."""

        phases = np.asarray(phases)

        # Now assign phases:
        rng = np.random.default_rng(seed=random_seed)
        phase_phase_type = rng.choice(
            a=self.num_phase_types,
            size=phases.size,
            p=self.target_phase_type_fractions,
        )
        for type_idx, phase_type in enumerate(self.phase_types):

            phase_idx_i = np.where(phase_phase_type == type_idx)[0]

            if phase_type.orientations is not None:
                num_oris_i = phase_type.orientations.shape[0]
                num_phases_i = len(phase_idx_i)
                if num_oris_i < num_phases_i:
                    raise ValueError(
                        f"Insufficient number of orientations ({num_oris_i}) for phase type "
                        f"{type_idx} with {num_phases_i} phases."
                    )
                elif num_oris_i > num_phases_i:
                    # select a subset randomly:
                    oris_i_idx = rng.choice(a=num_oris_i, size=num_phases_i)
                    phase_type.orientations = phase_type.orientations[oris_i_idx]

            phase_type.phases = phases[phase_idx_i]


class PhaseTypeDefinition:
    """Class to represent a type of phase (i.e. grain) within a CIPHER material.

    Attributes
    ----------
    material : MaterialDefinition
        Material to which this phase type belongs.
    type_label : str
        To distinguish between multiple phase types that all belong to the same material.
    target_type_fraction : float
    phases : ndarray of shape (N,) of int
        Phases that belong to this phase type.
    orientations : ndarray of shape (N, 4) of float
        Quaternion orientations for each phase.
    """

    def __init__(
        self,
        type_label=None,
        target_type_fraction=None,
        phases=None,
        orientations=None,
    ):
        self.type_label = type_label
        self.target_type_fraction = target_type_fraction
        self.phases = np.asarray(phases) if phases is not None else phases
        self.orientations = orientations

        self._material = None

        if self.phases is not None and self.target_type_fraction is not None:
            raise ValueError("Cannot specify both `phases` and `target_type_fraction`.")

        if orientations is not None and phases is None:
            raise ValueError("If specifying `orientations`, must also specify `phases`.")

    @property
    def material(self):
        return self._material

    @property
    def name(self):
        return self.material.name + (f"-{self.type_label}" if self.type_label else "")

    def to_JSON(self, keep_arrays=False):
        data = {
            "type_label": self.type_label,
            "phases": self.phases,
            "orientations": self.orientations,
        }
        if not keep_arrays:
            data["phases"] = data["phases"].tolist()
            if self.orientations is not None:
                data["orientations"] = data["orientations"].tolist()

        return data

    @classmethod
    def from_JSON(cls, data):
        data = {
            "type_label": data["type_label"],
            "phases": np.array(data["phases"]),
            "orientations": np.array(data["orientations"])
            if data["orientations"] is not None
            else None,
        }
        return cls(**data)


@dataclass
class CIPHERInput:
    geometry: CIPHERGeometry
    components: List
    outputs: List
    solution_parameters: Dict

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
    def from_JSON(cls, data):
        data = {
            "geometry": CIPHERGeometry.from_JSON(data["geometry"]),
            "components": data["components"],
            "outputs": data["outputs"],
            "solution_parameters": data["solution_parameters"],
        }
        return cls(**data)

    @classmethod
    def from_input_YAML_file(cls, path):
        """Generate a CIPHERInput object from a CIPHER input YAML file."""

        with Path(path).open("rt") as fp:
            file_str = "".join(fp.readlines())

        return cls.from_input_YAML_str(file_str)

    @classmethod
    def read_input_YAML_file(cls, path):

        with Path(path).open("rt") as fp:
            file_str = "".join(fp.readlines())

        return cls.read_input_YAML_string(file_str)

    @staticmethod
    def read_input_YAML_string(file_str, parse_interface_map=True):

        yaml = YAML(typ="safe")
        data = yaml.load(file_str)

        header = data["header"]
        grid_size = header["grid"]
        size = header["size"]
        num_phases = header["n_phases"]

        voxel_phase = decompress_1D_array_string(data["mappings"]["voxel_phase_mapping"])
        voxel_phase = voxel_phase.reshape(grid_size, order="F") - 1

        unique_phase_IDs = np.unique(voxel_phase)
        assert len(unique_phase_IDs) == num_phases

        interface_map = None
        if parse_interface_map:
            interface_map = decompress_1D_array_string(
                data["mappings"]["interface_mapping"]
            )
            interface_map = interface_map.reshape((num_phases, num_phases)) - 1
            interface_map[np.tril_indices(num_phases)] = -1  # only need one half

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
    def from_input_YAML_str(cls, file_str):
        """Generate a CIPHERInput object from a CIPHER input YAML file string."""

        yaml_dat = cls.read_input_YAML_string(file_str)
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
        )

        attrs = {
            "geometry": geom,
            "components": yaml_dat["header"]["components"],
            "outputs": yaml_dat["header"]["outputs"],
            "solution_parameters": dict(yaml_dat["solution_parameters"]),
        }

        return cls(**attrs)

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

    def write_yaml(self, path):
        """Write the CIPHER input YAML file."""

        self.geometry._validate_interface_map()

        cipher_input_data = {
            "header": self.get_header(),
            "solution_parameters": dict(sorted(self.solution_parameters.items())),
            "material": {
                k: copy.deepcopy(v) for k, v in self.material_properties.items()
            },
            "interface": {k: copy.deepcopy(v) for k, v in self.get_interfaces().items()},
            "mappings": {
                "phase_material_mapping": LiteralScalarString(
                    compress_1D_array_string(self.geometry.phase_material + 1) + "\n"
                ),
                "voxel_phase_mapping": LiteralScalarString(
                    compress_1D_array_string(
                        self.geometry.voxel_phase.flatten(order="F") + 1
                    )
                    + "\n"
                ),
                "interface_mapping": LiteralScalarString(
                    compress_1D_array_string(
                        self.geometry.interface_map_int.flatten() + 1
                    )
                    + "\n"
                ),
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
        energy_range,
        mobility_range,
        theta_max,
        n=4,
        B=5,
        bin_width=5,
        degrees=True,
    ):
        base_defn, phase_pairs = self.geometry.remove_interface(base_interface_name)
        if self.geometry.misorientation_matrix is None:
            misori_matrix = self.geometry.get_misorientation_matrix()
        else:
            misori_matrix = self.geometry.misorientation_matrix

        print(f"{bin_width=}")

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
        num_bins = misori_bins.size
        energy_bins = np.linspace(*energy_range, num=num_bins)
        mobility_bins = np.linspace(*mobility_range, num=num_bins)

        theta = (misori_bins + (bin_width / 2))[:-1]
        print(f"{misori_bins=}")
        print(f"{theta=}")

        energy = (
            read_shockley(
                theta=theta,
                E_max=(energy_range[1] - energy_range[0]),
                theta_max=theta_max,
                degrees=degrees,
            )
            + energy_range[0]
        )
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

                new_e0 = energy[bin_idx_i - 1].item()
                set_by_path(root=props, path=("energy", "e0"), value=new_e0)

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
