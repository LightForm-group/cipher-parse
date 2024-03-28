from damask import Orientation
import pyvista as pv
import numpy as np
import plotly.express as px

from cipher_parse.material import MaterialDefinition
from cipher_parse.interface import InterfaceDefinition
from cipher_parse.discrete_voronoi import DiscreteVoronoi
from cipher_parse.voxel_map import VoxelMap
from cipher_parse.errors import (
    GeometryDuplicateMaterialNameError,
    GeometryExcessTargetVolumeFractionError,
    GeometryMissingPhaseAssignmentError,
    GeometryNonUnitTargetVolumeFractionError,
    GeometryUnassignedPhasePairInterfaceError,
    GeometryVoxelPhaseError,
)
from cipher_parse.utilities import generate_interface_energies_plot
from cipher_parse.quats import (
    compute_misorientation_matrix,
    quat_angle_between,
    compute_misorientation_matrix_damask,
)


class CIPHERGeometry:
    def __init__(
        self,
        materials,
        interfaces,
        size,
        seeds=None,
        voxel_phase=None,
        voxel_map=None,
        is_periodic=False,
        random_seed=None,
        allow_missing_phases=False,
        quiet=False,
        time=None,
        increment=None,
        incremental_data_idx=None,
    ):
        """
        Parameters
        ----------

        allow_missing_phases : bool, optional
            If True, allow references to phases that do not appear in the voxel map.

        """

        if sum(i is not None for i in (voxel_phase, voxel_map)) != 1:
            raise ValueError(f"Specify exactly one of `voxel_phase` and `voxel_map`")
        if voxel_map is None:
            voxel_map = VoxelMap(
                region_ID=voxel_phase,
                size=size,
                is_periodic=is_periodic,
                quiet=quiet,
            )
        else:
            voxel_phase = voxel_map.region_ID
            is_periodic = voxel_map.is_periodic

        self._interfaces = None
        self._is_periodic = is_periodic

        self.voxel_map = voxel_map
        self.voxel_phase = voxel_phase
        self.seeds = np.asarray(seeds)
        self.materials = materials
        self.interfaces = interfaces
        self.random_seed = random_seed
        self.size = np.asarray(size)
        self.allow_missing_phases = allow_missing_phases
        self.time = time
        self.increment = increment
        self.incremental_data_idx = incremental_data_idx

        for i in self.materials:
            i._geometry = self

        for idx, i in enumerate(self.interfaces):
            i.index = idx

        if self.size.size != self.dimension:
            raise ValueError(
                f"`size` ({self.size}) implies {self.size.size} dimensions, but "
                f"`voxel_phase` implies {self.voxel_phase.ndim} dimensions."
            )

        all_phases = self.present_phases
        self._num_phases = all_phases.size

        if not allow_missing_phases:
            if not np.all(all_phases == np.arange(self.num_phases)):
                raise GeometryVoxelPhaseError(
                    "`voxel_phase` must be an array of consecutive integers starting from "
                    "zero."
                )

        if len(set(self.material_names)) < self.num_materials:
            raise GeometryDuplicateMaterialNameError(
                f"Repeated material names exist in the materials definitions: "
                f"{self.material_names!r}."
            )

        self._ensure_phase_assignment(random_seed)

        self._phase_material = self._get_phase_material()
        self._validate_interfaces()
        self._check_interface_phase_pairs()

        self._phase_phase_type = self._get_phase_phase_type()

        # assigned by calculate_* methods on first call to corresponding get_* methods:
        self._phase_voxels = None
        self._phase_num_voxels = None
        self._phase_voxel_indices = None
        self._phase_voxel_coordinates = None
        self._phase_voxel_centroids = None
        self._grain_boundaries = None
        self._grain_boundary_centroids = None

        self._interface_map = self._get_interface_map(quiet=quiet)
        self._validate_interface_map()  # TODO: add setter to interface map

        self.phase_orientation = self._get_phase_orientation()

        # assigned by `get_misorientation_matrix`:
        self._misorientation_matrix = None
        self._misorientation_matrix_is_degrees = None

    def __eq__(self, other):
        # Note we don't check seeds (not stored in YAML file)
        if not isinstance(other, self.__class__):
            return False
        if (
            self.materials == other.materials
            and self.interfaces == other.interfaces
            and np.all(self.size == self.size)
            and np.all(self.random_seed == self.random_seed)
            and np.all(self.voxel_phase == self.voxel_phase)
        ):
            return True
        return False

    def _validate_interfaces(self):
        int_names = self.interface_names
        if len(set(int_names)) < len(int_names):
            raise ValueError(
                f"Multiple interfaces have the same name (i.e. "
                f"phase-type-pair and type-label combination)!"
            )

    def to_JSON(self, keep_arrays=False):
        data = {
            "materials": [i.to_JSON(keep_arrays) for i in self.materials],
            "interfaces": [i.to_JSON(keep_arrays) for i in self.interfaces],
            "size": self.size,
            "seeds": self.seeds,
            "voxel_phase": self.voxel_phase,
            "is_periodic": self.is_periodic,
            "random_seed": self.random_seed,
            "misorientation_matrix": self.misorientation_matrix,
            "misorientation_matrix_is_degrees": self.misorientation_matrix_is_degrees,
            "allow_missing_phases": self.allow_missing_phases,
            "grain_boundaries": self._grain_boundaries,
            "time": self.time,
            "increment": self.increment,
            "incremental_data_idx": self.incremental_data_idx,
        }
        if not keep_arrays:
            data["size"] = data["size"].tolist()
            data["seeds"] = data["seeds"].tolist()
            data["voxel_phase"] = data["voxel_phase"].tolist()
            if data["misorientation_matrix"] is not None:
                data["misorientation_matrix"] = data["misorientation_matrix"].tolist()
            for phase_pair in data.get("grain_boundaries") or []:
                GB = data["grain_boundaries"][phase_pair]
                data["grain_boundaries"][phase_pair] = {
                    "interface_idx": GB["interface_idx"],
                    "voxel_indices": list(i.tolist() for i in GB["voxel_indices"]),
                    "voxel_coordinates": GB["voxel_coordinates"].tolist(),
                    "centroid": GB["centroid"].tolist(),
                }

        return data

    @classmethod
    def from_JSON(cls, data, quiet=True):
        data_init = {
            "materials": [MaterialDefinition.from_JSON(i) for i in data["materials"]],
            "interfaces": [InterfaceDefinition.from_JSON(i) for i in data["interfaces"]],
            "size": np.array(data["size"]),
            "seeds": np.array(data["seeds"]),
            "voxel_phase": np.array(data["voxel_phase"]),
            "is_periodic": data.get("is_periodic", False),
            "random_seed": data["random_seed"],
            "allow_missing_phases": data.get("allow_missing_phases", False),
            "time": data.get("time"),
            "increment": data.get("increment"),
            "incremental_data_idx": data.get("incremental_data_idx"),
        }
        GBs = {}
        for phase_pair in data.get("grain_boundaries") or []:
            GB = data["grain_boundaries"][phase_pair]
            GBs[phase_pair] = {
                "interface_idx": GB["interface_idx"],
                "voxel_indices": tuple(np.array(i) for i in GB["voxel_indices"]),
                "voxel_coordinates": np.array(GB["voxel_coordinates"]),
                "centroid": np.array(GB["centroid"]),
            }
        obj = cls(**data_init, quiet=quiet)
        if data["misorientation_matrix"] is not None:
            obj._misorientation_matrix = np.array(data["misorientation_matrix"])
        obj._grain_boundaries = GBs or None
        return obj

    @property
    def present_phases(self):
        return np.unique(self.voxel_phase)

    @property
    def missing_phases(self):
        return np.array(list(set(self.known_phases) - set(self.present_phases)))

    @property
    def known_phases(self):
        phases = []
        for mat in self.materials:
            phases.append(mat.phases)
        return np.concatenate(phases)

    @property
    def interfaces(self):
        return self._interfaces

    @property
    def is_periodic(self):
        return self._is_periodic

    @interfaces.setter
    def interfaces(self, interfaces):
        self._interfaces = interfaces
        self._validate_interfaces()

    @property
    def misorientation_matrix(self):
        return self._misorientation_matrix

    @property
    def misorientation_matrix_is_degrees(self):
        return self._misorientation_matrix_is_degrees

    def get_phase_voxels(self):
        if self._phase_voxels is None:
            self._calculate_phase_voxels()
        return self._phase_voxels

    def get_phase_num_voxels(self):
        if self._phase_num_voxels is None:
            self._calculate_phase_num_voxels()
        return self._phase_num_voxels

    def get_phase_voxel_indices(self):
        if self._phase_voxel_indices is None:
            self._calculate_phase_voxel_indices()
        return self._phase_voxel_indices

    def get_phase_voxel_coordinates(self):
        if self._phase_voxel_coordinates is None:
            self._calculate_phase_voxel_coordinates()
        return self._phase_voxel_coordinates

    def get_phase_voxel_centroids(self):
        if self._phase_voxel_centroids is None:
            self._calculate_phase_voxel_centroids()
        return self._phase_voxel_centroids

    def get_grain_boundaries(self):
        if self._grain_boundaries is None:
            self._calculate_grain_boundaries()
        return self._grain_boundaries

    def get_grain_boundary_centroids(self):
        if self._grain_boundary_centroids is None:
            self._calculate_grain_boundary_centroids()
        return self._grain_boundary_centroids

    def _calculate_phase_voxels(self):
        self._phase_voxels = [
            self.voxel_phase == phase_idx for phase_idx in range(self.num_known_phases)
        ]

    def _calculate_phase_num_voxels(self):
        self._phase_num_voxels = np.array([np.sum(i) for i in self.get_phase_voxels()])

    def _calculate_phase_voxel_indices(self):
        self._phase_voxel_indices = [np.where(i) for i in self.get_phase_voxels()]

    def _calculate_phase_voxel_coordinates(self):
        self._phase_voxel_coordinates = [
            self.voxel_map.coordinates[i] for i in self.get_phase_voxel_indices()
        ]

    def _calculate_phase_voxel_centroids(self):
        self._phase_voxel_centroids = np.array(
            [
                np.mean(i, axis=0) if i.size else np.ones((self.dimension,)) * np.nan
                for i in self.get_phase_voxel_coordinates()
            ]
        )

    def _calculate_grain_boundaries(self):
        grain_boundaries = {}
        tot_num_calcs = sum(i.phase_pairs.shape[0] for i in self.interfaces)
        calc_count = 0
        report_each_pc = 5
        num_iter_per_report = np.ceil(tot_num_calcs * report_each_pc / 100)
        print(f"Identifying grain boundaries...", flush=True)
        for int_idx, interface in enumerate(self.interfaces):
            for phase_pair in interface.phase_pairs:
                calc_count += 1
                if calc_count % num_iter_per_report == 0:
                    frac_done = (1 + calc_count) / tot_num_calcs * 100
                    print(f"Identifying grain boundaries: {frac_done:.0f}%.", flush=True)
                is_GB = np.any(np.all(self.neighbour_list == phase_pair[:, None], axis=0))
                if is_GB:
                    vox_bool = self.voxel_map.get_region_boundary_voxels(
                        phase_pair[0], phase_pair[1]
                    )
                    vox_idx = np.where(vox_bool)

                    # remove edge GB voxels:
                    for gs_idx, i in enumerate(self.grid_size):
                        not_edge_idx = np.logical_not(
                            np.logical_or(vox_idx[gs_idx] == 0, vox_idx[gs_idx] == i - 1)
                        )
                        vox_idx = list(vox_idx)
                        for j_idx, _ in enumerate(vox_idx):
                            vox_idx[j_idx] = vox_idx[j_idx][not_edge_idx]

                    vox_idx = tuple(vox_idx)
                    if not vox_idx[0].size:
                        continue

                    vox_coords = self.voxel_map.coordinates[vox_idx]
                    GB_centroid = np.mean(vox_coords, axis=0)

                    grain_boundaries[(phase_pair[0], phase_pair[1])] = {
                        "interface_idx": int_idx,
                        "voxel_indices": vox_idx,
                        "voxel_coordinates": vox_coords,
                        "centroid": GB_centroid,
                    }
        print(f"Finished grain boundaries.", flush=True)
        self._grain_boundaries = grain_boundaries

    def _calculate_grain_boundary_centroids(self):
        self._grain_boundary_centroids = np.concatenate(
            [i["centroid"][None] for i in self.get_grain_boundaries().values()], axis=0
        )

    def _ensure_phase_assignment(self, random_seed):
        is_mat_phases = [i.phases is not None for i in self.materials]
        is_mat_vol_frac = [i is not None for i in self.target_material_volume_fractions]
        is_mixed = any(is_mat_phases) and any(is_mat_vol_frac)

        if is_mixed or (any(is_mat_phases) and not all(is_mat_phases)):
            raise GeometryMissingPhaseAssignmentError(
                f"Specify either: all phases explicitly (via the material definition "
                f"`phases`, or the constituent phase type definition `phases`), or "
                f"specify zero or more target volume fractions for the material "
                f"definitions."
            )

        if not any(is_mat_phases):
            self._assign_phases_by_volume_fractions(is_mat_vol_frac, random_seed)

    def _check_interface_phase_pairs(self):
        """If interfaces have phase-pairs specified, check these are consistent with
        the specified phases of associated material."""

        for i in self.interfaces:
            # assign materials as well as phase_types if materials not assigned to
            # interface:
            if not i.materials:
                i_mats = []
                # find which material each referenced phase type belongs to:
                for j in i.phase_types:
                    mat_j = [k.material.name for k in self.phase_types if k.name == j][0]
                    i_mats.append(mat_j)
                i.materials = tuple(i_mats)

            if i.phase_pairs.size:
                mats_idx = np.sort([self.material_names.index(j) for j in i.materials])
                phase_pairs_material = self.phase_material[i.phase_pairs]
                phase_pairs_mat_srt = np.sort(phase_pairs_material, axis=1)
                if not np.all(np.all(phase_pairs_mat_srt == mats_idx, axis=1)):
                    raise ValueError(
                        f"Phase pairs specified for interface {i.name!r} are not "
                        f"consistent with phases specified for the interface materials "
                        f"{i.materials[0]!r} and {i.materials[1]!r}."
                    )  # TODO: test raise

    def _assign_phases_by_volume_fractions(self, is_mat_vol_frac, random_seed):
        # Assign via target volume fractions.
        num_unassigned_vol = self.num_materials - sum(is_mat_vol_frac)
        assigned_vol = sum(i or 0.0 for i in self.target_material_volume_fractions)
        if num_unassigned_vol:
            frac = (1.0 - assigned_vol) / num_unassigned_vol
            if frac <= 0.0:
                raise GeometryExcessTargetVolumeFractionError(
                    f"All material target volume fractions must sum to one, but "
                    f"assigned target volume fractions sum to {assigned_vol} with "
                    f"{num_unassigned_vol} outstanding unassigned material volume "
                    f"fraction(s)."
                )
        for i in self.materials:
            if i.target_volume_fraction is None:
                i.target_volume_fraction = frac

        assigned_vol = sum(self.target_material_volume_fractions)
        if not np.isclose(assigned_vol, 1.0):
            raise GeometryNonUnitTargetVolumeFractionError(
                f"All material target volume fractions must sum to one, but target "
                f"volume fractions sum to {assigned_vol}."
            )

        # Now assign phases:
        rng = np.random.default_rng(seed=random_seed)
        phase_material = rng.choice(
            a=self.num_materials,
            size=self.num_phases,
            p=self.target_material_volume_fractions,
        )
        for mat_idx, mat in enumerate(self.materials):
            mat_phases = np.where(phase_material == mat_idx)[0]
            mat.assign_phases(mat_phases)

    def _get_phase_material(self):
        phase_material = np.ones(self.num_known_phases) * np.nan
        all_phase_idx = []
        for mat_idx, mat in enumerate(self.materials):
            try:
                phase_material[mat.phases] = mat_idx
                all_phase_idx.append(mat.phases)
            except IndexError:
                raise ValueError(
                    f"Material {mat.name!r} phases indices {mat.phases} are invalid, "
                    f"given the number of phases ({self.num_phases})."
                )
        if np.any(np.isnan(phase_material)):
            raise ValueError(
                "Not all phases are accounted for in the phase type definitions."
            )  # TODO: test raise

        # check all phase indices form a consequtive range:
        num_phases_range = set(np.arange(self.num_phases))
        known_phases = set(np.hstack(all_phase_idx))
        miss_phase_idx = num_phases_range - known_phases
        bad_phase_idx = known_phases - num_phases_range
        if miss_phase_idx:
            raise ValueError(
                f"Missing phase indices: {miss_phase_idx}. Bad phase indices: "
                f"{bad_phase_idx}"
            )  # TODO: test raise

        return phase_material.astype(int)

    def _get_phase_phase_type(self):
        phase_phase_type = np.ones(self.num_known_phases) * np.nan
        for phase_type_idx, phase_type in enumerate(self.phase_types):
            phase_phase_type[phase_type.phases] = phase_type_idx
        if np.any(np.isnan(phase_phase_type)):
            raise RuntimeError("Not all phases accounted for!")  # TODO: test raise?
        return phase_phase_type.astype(int)

    def _get_phase_orientation(self):
        """Get the orientation of each phase, if specified in the phase-type."""
        phase_ori = np.ones((self.num_known_phases, 4), dtype=float) * np.nan
        for phase_type in self.phase_types:
            for type_idx, phase_i in enumerate(phase_type.phases):
                if phase_type.orientations is not None:
                    phase_ori[phase_i] = phase_type.orientations[type_idx]
        return phase_ori

    def get_interface_map_indices(self, phase_type_A, phase_type_B):
        """Get an array of integer indices that index the (upper triangle of the) 2D
        symmetric interface map array, corresponding to a given material pair."""

        # First get phase indices belonging to the two phase types:
        ptypes = {i.name: i for i in self.phase_types}
        ptA_phases = ptypes[phase_type_A].phases
        ptB_phases = ptypes[phase_type_B].phases

        A_idx = np.repeat(ptA_phases, ptB_phases.shape[0])
        B_idx = np.tile(ptB_phases, ptA_phases.shape[0])

        map_idx = np.vstack((A_idx, B_idx))
        map_idx_srt = np.sort(map_idx, axis=0)  # map onto upper triangle
        map_idx_uniq = np.unique(map_idx_srt, axis=1)  # get unique pairs only

        # remove diagonal elements (a phase can't have an interface with itself)
        map_idx_non_trivial = map_idx_uniq[:, map_idx_uniq[0] != map_idx_uniq[1]]

        return map_idx_non_trivial

    def _get_interface_map(self, upper_tri_only=False, quiet=False):
        """Generate the num_phases by num_phases symmetric matrix that maps each phase-pair
        to an interface index."""

        if not quiet:
            print("Finding interface map matrix...", end="")

        int_map = (
            np.ones((self.num_known_phases, self.num_known_phases), dtype=int) * np.nan
        )

        ints_by_phase_type_pair = {}
        for int_def in self.interfaces:
            if int_def.phase_types not in ints_by_phase_type_pair:
                ints_by_phase_type_pair[int_def.phase_types] = []
            ints_by_phase_type_pair[int_def.phase_types].append(int_def)

        for pt_pair, int_defs in ints_by_phase_type_pair.items():
            names = [i.name for i in int_defs]
            if len(set(names)) < len(names):
                raise ValueError(
                    f"Multiple interface definitions for phase-type pair "
                    f"{pt_pair} have the same `type_label`."
                )
            type_fracs = [i.type_fraction for i in int_defs]
            any_frac_set = any(i is not None for i in type_fracs)
            manual_set = [i.is_phase_pairs_set for i in int_defs]
            any_manual_set = any(manual_set)
            all_manual_set = all(manual_set)
            if any_frac_set:
                if any_manual_set:
                    raise ValueError(
                        f"For interface {pt_pair}, specify phase pairs manually for all "
                        f"defined interfaces using `phase_pairs`, or specify `type_fraction`"
                        f"for all defined interfaces. You cannot mix them."
                    )

            all_phase_pairs = self.get_interface_map_indices(*pt_pair).T
            if any_manual_set:
                if not all_manual_set:
                    raise ValueError(
                        f"For interface {pt_pair}, specify phase pairs manually for all "
                        f"defined interfaces using `phase_pairs`, or specify `type_fraction`"
                        f"for all defined interfaces. You cannot mix them."
                    )

                # check that given phase_pairs combine to the set of all phase_pairs
                # for this material-material pair:
                all_given_phase_pairs = np.vstack([i.phase_pairs for i in int_defs])

                # sort by first-phase, then second-phase, for comparison:
                srt = np.lexsort(all_given_phase_pairs.T[::-1])
                all_given_phase_pairs = all_given_phase_pairs[srt]

                if all_given_phase_pairs.shape != all_phase_pairs.shape or not np.all(
                    all_given_phase_pairs == all_phase_pairs
                ):
                    raise ValueError(
                        f"Missing `phase_pairs` for interface {pt_pair}. The following "
                        f"phase pairs must all be included for this interface: "
                        f"{all_phase_pairs}"
                    )

                for int_i in int_defs:
                    phase_pairs_i = int_i.phase_pairs.T
                    if phase_pairs_i.size:
                        int_map[phase_pairs_i[0], phase_pairs_i[1]] = int_i.index

                        if not upper_tri_only:
                            int_map[phase_pairs_i[1], phase_pairs_i[0]] = int_i.index

            else:
                # set default type fractions if missing
                remainder_frac = 1 - sum(i for i in type_fracs if i is not None)
                if remainder_frac > 0:
                    num_missing_type_frac = sum(1 for i in type_fracs if i is None)
                    if num_missing_type_frac == 0:
                        raise ValueError(
                            f"For interface {pt_pair}, `type_fraction` for all "
                            f"defined interfaces must sum to one."
                        )
                    remainder_frac_each = remainder_frac / num_missing_type_frac
                    for i in int_defs:
                        if i.type_fraction is None:
                            i.type_fraction = remainder_frac_each

                type_fracs = [i.type_fraction for i in int_defs]
                if sum(type_fracs) != 1:
                    raise ValueError(
                        f"For interface {pt_pair}, `type_fraction` for all "
                        f"defined interfaces must sum to one."
                    )

                # assign phase_pairs according to type fractions:
                num_pairs = all_phase_pairs.shape[0]
                type_nums_each = [round(i * num_pairs) for i in type_fracs]
                type_nums = np.cumsum(type_nums_each)
                if num_pairs % 2 == 1:
                    type_nums += 1

                shuffle_idx = np.random.choice(num_pairs, size=num_pairs, replace=False)
                phase_pairs_shuffled = all_phase_pairs[shuffle_idx]
                phase_pairs_split = np.split(phase_pairs_shuffled, type_nums, axis=0)[:-1]
                for idx, int_i in enumerate(int_defs):
                    phase_pairs_i = phase_pairs_split[idx]
                    int_map[phase_pairs_i[:, 0], phase_pairs_i[:, 1]] = int_i.index
                    if not upper_tri_only:
                        int_map[phase_pairs_i[:, 1], phase_pairs_i[:, 0]] = int_i.index
                    int_i.phase_pairs = phase_pairs_i
                    int_i.type_fraction = None

        if not quiet:
            print("done!")

        return int_map

    @property
    def interface_map_int(self):
        """Get the interface map as an integer matrix, where NaNs are replaced by -2."""
        int_map = np.copy(self.interface_map)
        int_map[np.isnan(int_map)] = -2
        return int_map.astype(int)

    def get_interface_idx(self):
        """Get the interface index associated with each voxel."""
        return self.voxel_map.get_interface_idx(self.interface_map_int)

    def get_interface_misorientation(self):
        return self.voxel_map.get_interface_idx(self.misorientation_matrix)

    def _modify_interface_map(self, phase_A, phase_B, interface_idx):
        """
        Parameters
        ----------
        phase_A : ndarray
        phase_B : ndarray
        interface_idx : ndarray

        """
        if interface_idx not in range(len(self.interfaces)):
            raise ValueError(f"Interface index {interface_idx} invalid.")
        self._interface_map[phase_A, phase_B] = interface_idx
        self._interface_map[phase_B, phase_A] = interface_idx

    def _validate_interface_map(self):
        # check no missing interfaces:
        int_map_indices = np.triu_indices_from(self.interface_map, k=1)
        int_is_nan = np.isnan(self.interface_map[int_map_indices])
        phase_idx_int_is_nan = np.vstack(int_map_indices)[:, int_is_nan]
        if phase_idx_int_is_nan.size:
            raise GeometryUnassignedPhasePairInterfaceError(
                f"The following phase-pairs have not been assigned an interface "
                f"definition: {phase_idx_int_is_nan}."
            )

    def get_misorientation_matrix(self, degrees=True, overwrite=False, method="defdap"):
        """Given phase type definitions that include orientation lists, get the
        misorientation matrix between all pairs."""

        if self.misorientation_matrix is not None and not overwrite:
            print(
                "Misorientation matrix is already set. Use `overwrite=True` to recompute."
            )
            return

        all_oris = np.ones((self.num_phases, 4)) * np.nan
        for i in self.phase_types:
            all_oris[i.phases] = i.orientations

        if np.any(np.isnan(all_oris)):
            raise RuntimeError(
                "Not all orientations are accounted for in the phase type definitions."
            )

        if method == "defdap":
            misori_matrix = compute_misorientation_matrix(all_oris, degrees)
        elif method == "damask":
            misori_matrix = compute_misorientation_matrix_damask(all_oris, degrees)

        self._misorientation_matrix = misori_matrix
        self._misorientation_matrix_is_degrees = degrees

        return misori_matrix

    def get_pyvista_grid(self):
        """Experimental!"""

        grid = pv.UniformGrid()

        grid.dimensions = self.grid_size_3D + 1  # +1 to inject values on cell data
        grid.spacing = self.size_3D / self.grid_size_3D
        return grid

    @staticmethod
    def get_unique_random_seeds(num_phases, size, grid_size, random_seed=None):
        return DiscreteVoronoi.get_unique_random_seeds(
            num_regions=num_phases,
            size=size,
            grid_size=grid_size,
            random_seed=random_seed,
        )

    @staticmethod
    def assign_phase_material_randomly(
        num_materials,
        num_phases,
        volume_fractions,
        random_seed=None,
    ):
        print(
            "Randomly assigning phases to materials according to volume_fractions...",
            end="",
        )
        rng = np.random.default_rng(seed=random_seed)
        phase_material = rng.choice(
            a=num_materials,
            size=num_phases,
            p=volume_fractions,
        )
        print("done!")
        return phase_material

    @classmethod
    def from_voronoi(
        cls,
        interfaces,
        materials,
        grid_size,
        size,
        seeds=None,
        num_phases=None,
        random_seed=None,
        is_periodic=False,
    ):
        if sum(i is not None for i in (seeds, num_phases)) != 1:
            raise ValueError(f"Specify exactly one of `seeds` and `num_phases`")

        if seeds is None:
            vor_map = DiscreteVoronoi.from_random(
                num_regions=num_phases,
                grid_size=grid_size,
                size=size,
                is_periodic=is_periodic,
                random_seed=random_seed,
            )
            seeds = vor_map.seeds

        else:
            vor_map = DiscreteVoronoi.from_seeds(
                region_seeds=seeds,
                grid_size=grid_size,
                size=size,
                is_periodic=is_periodic,
            )

        return cls(
            voxel_map=vor_map,
            materials=materials,
            interfaces=interfaces,
            size=size,
            seeds=seeds,
            random_seed=random_seed,
        )

    @classmethod
    def from_seed_voronoi(
        cls,
        seeds,
        interfaces,
        materials,
        grid_size,
        size,
        random_seed=None,
        is_periodic=False,
    ):
        return cls.from_voronoi(
            interfaces=interfaces,
            materials=materials,
            grid_size=grid_size,
            size=size,
            seeds=seeds,
            random_seed=random_seed,
            is_periodic=is_periodic,
        )

    @classmethod
    def from_random_voronoi(
        cls,
        num_phases,
        interfaces,
        materials,
        grid_size,
        size,
        random_seed=None,
        is_periodic=False,
    ):
        return cls.from_voronoi(
            interfaces=interfaces,
            materials=materials,
            grid_size=grid_size,
            size=size,
            num_phases=num_phases,
            random_seed=random_seed,
            is_periodic=is_periodic,
        )

    @property
    def voxel_phase_3D(self):
        if self.dimension == 3:
            return self.voxel_phase
        else:
            return self.voxel_phase.T[:, :, None]

    @property
    def voxel_material_3D(self):
        if self.dimension == 3:
            return self.voxel_material
        else:
            return self.voxel_material.T[:, :, None]

    @property
    def voxel_interface_idx(self):
        return self.get_interface_idx()

    @property
    def voxel_interface_idx_3D(self):
        int_idx = self.voxel_interface_idx
        if self.dimension == 3:
            return int_idx
        else:
            return int_idx.T[:, :, None]

    @property
    def voxel_interface(self):
        return self.voxel_map.get_interface_voxels()

    @property
    def voxel_interface_3D(self):
        int_vox = self.voxel_interface
        if self.dimension == 3:
            return int_vox
        else:
            return int_vox.T[:, :, None]

    @property
    def voxel_phase_neighbours_3D(self):
        if self.dimension == 3:
            return self.voxel_phase_neighbours
        else:
            return self.voxel_phase_neighbours.T[:, :, None]

    def get_slice(
        self,
        slice_index=0,
        normal_dir="z",
        data_label="phase",
        include=None,
        misorientation_matrix=None,
    ):
        allowed_data = [
            "phase",
            "material",
            "interface",
            "interface_idx",
            "phase_neighbours",
            "grain_boundaries",
            "GB_misorientation",
            "IPF_z",
        ]
        if data_label not in allowed_data:
            raise ValueError(f"`data_label` must be one of: {allowed_data}.")

        if data_label == "phase":
            data = self.voxel_phase_3D
        elif data_label == "material":
            data = self.voxel_material_3D
        elif data_label == "interface":
            data = self.voxel_interface_3D
        elif data_label == "interface_idx":
            data = self.voxel_interface_idx_3D
        elif data_label == "phase_neighbours":
            data = self.voxel_phase_neighbours_3D
        elif data_label == "grain_boundaries":
            data = self.get_grain_boundary_map(as_3D=True)
        elif data_label == "GB_misorientation":
            data = self.voxel_map.get_interface_idx(misorientation_matrix, as_3D=True)
        elif data_label == "IPF_z":
            data = self.get_voxel_IPF(IPF_dir=None, as_3D=True)

        data = np.copy(data)
        if normal_dir == "x":
            data = data[slice_index, :, :]
        elif normal_dir == "y":
            data = data[:, slice_index, :]
        elif normal_dir == "z":
            data = data[:, :, slice_index]

        if include:
            include_mask = data == include[0]
            for i in include[1:]:
                include_mask = np.logical_or(include_mask, data == i)

            data[~include_mask] = -10

        return data

    def show_slice(
        self,
        slice_index=0,
        normal_dir="z",
        data_label="phase",
        include=None,
        phase_centroids=False,
        discrete_colours=None,
        layout_args=None,
        **kwargs,
    ):
        if "misorientation_matrix" in kwargs:
            slice_dat = self.get_slice(
                slice_index,
                normal_dir,
                data_label,
                include,
                misorientation_matrix=kwargs.pop("misorientation_matrix"),
            )
        else:
            slice_dat = self.get_slice(
                slice_index,
                normal_dir,
                data_label,
                include,
            )

        if discrete_colours:
            slice_RGB = np.tile(slice_dat[..., None], (1, 1, 3)).astype(float)
            for k, v in discrete_colours.items():
                slice_RGB[np.where(np.all(slice_RGB == k, axis=2))] = v
            slice_dat = slice_RGB

        fig = px.imshow(
            slice_dat,
            color_continuous_scale="viridis",
            **kwargs,
        )
        if phase_centroids:
            # TODO fix for 3D?
            cents = self.get_phase_voxel_centroids()
            fig.add_scatter(
                x=cents[:, 1],
                y=cents[:, 0],
                text=np.arange(cents.shape[0]),
                mode="markers",
            )
        fig.update_layout(layout_args or {})
        return fig

    def show(self):
        """Experimental!"""

        print("WARNING: experimental!")

        grid = self.get_pyvista_grid()

        grid.cell_data["interface_idx"] = self.voxel_interface_idx_3D.flatten(order="F")
        grid.cell_data["material"] = self.voxel_material_3D.flatten(order="F")
        grid.cell_data["phase"] = self.voxel_phase_3D.flatten(order="F")

        # TODO: fix plotter to show multiple cell data
        pl = pv.Plotter(notebook=True)
        pl.add_mesh(grid)
        pl.show()

    def write_VTK(self, path):
        grid = self.get_pyvista_grid()

        grid.cell_data["interface_idx"] = self.voxel_interface_idx_3D.flatten(order="F")
        grid.cell_data["material"] = self.voxel_material_3D.flatten(order="F")
        grid.cell_data["phase"] = self.voxel_phase_3D.flatten(order="F")

        grid.save(path)

    @property
    def dimension(self):
        return self.voxel_map.dimension

    @property
    def grid_size(self):
        return self.voxel_map.grid_size

    @property
    def grid_size_3D(self):
        return self.voxel_map.grid_size_3D

    @property
    def size_3D(self):
        return self.voxel_map.size_3D

    @property
    def spacing(self):
        return self.voxel_map.spacing

    @property
    def spacing_3D(self):
        return self.voxel_map.spacing_3D

    def get_coordinates(self):
        return self.voxel_map.get_coordinates()

    @property
    def voxel_phase_neighbours(self):
        return self.voxel_map.neighbour_voxels

    @property
    def neighbour_list(self):
        return self.voxel_map.neighbour_list

    @property
    def interface_map(self):
        """Get the num_phases-by-num_phases matrix of interface indices."""
        return self._interface_map

    @property
    def interface_names(self):
        return [i.name for i in self.interfaces]

    @property
    def material_properties(self):
        return {mat.name: mat.properties for mat in self.materials}

    @property
    def num_voxels(self):
        return np.product(self.voxel_phase.size)  # TODO: change to voxel_map.num_voxels?

    @property
    def phase_voxels(self):
        return self._phase_voxels

    @property
    def phase_num_voxels(self):
        return self._phase_num_voxels

    @property
    def num_phases(self):
        return self._num_phases

    @property
    def num_known_phases(self):
        return len(self.known_phases)

    @property
    def num_materials(self):
        return len(self.materials)

    @property
    def material_names(self):
        return [i.name for i in self.materials]

    @property
    def target_material_volume_fractions(self):
        return [i.target_volume_fraction for i in self.materials]

    @property
    def phase_types(self):
        """Get all phase types across all materials."""
        return [j for i in self.materials for j in i.phase_types]

    @property
    def phase_material(self):
        """Get the material index of each phase."""
        return self._phase_material

    @property
    def phase_phase_type(self):
        """Get the phase type index of each phase."""
        return self._phase_phase_type

    @property
    def phase_orientation(self):
        """Get the orientation quaternion of each phase."""
        return self._phase_orientation

    @phase_orientation.setter
    def phase_orientation(self, oris):
        """Set the orientation quaternion of each phase, via material phase types."""

        shape = (self.num_known_phases, 4)
        if oris.shape != shape:
            raise ValueError(f"Phase orientations must have shape {shape!r}")

        max_idx = 0
        for phase_type in self.phase_types:
            ori_idx = np.arange(max_idx, max_idx + phase_type.num_phases)
            max_idx = ori_idx[-1] + 1
            phase_type.orientations = oris[ori_idx]

        self._phase_orientation = self._get_phase_orientation()

    @property
    def voxel_material(self):
        """Get the material index of each voxel."""
        return self.phase_material[self.voxel_phase]

    @property
    def voxel_phase_type(self):
        """Get the phase type index of each voxel."""
        return self.phase_phase_type[self.voxel_phase]

    @property
    def voxel_orientation(self):
        """Get the quaternion of each voxel."""
        return self.phase_orientation[self.voxel_phase]

    @property
    def material_num_voxels(self):
        mat_num_vox = []
        for i in self.materials:
            num_vox_i = sum(self.get_phase_num_voxels()[j] for j in i.phases)
            mat_num_vox.append(num_vox_i)
        return np.array(mat_num_vox)

    @property
    def phase_type_num_voxels(self):
        return np.array(
            [np.sum(self.get_phase_num_voxels()[i.phases]) for i in self.phase_types]
        )

    @property
    def material_volume_fractions(self):
        return np.array([i / self.num_voxels for i in self.material_num_voxels])

    @property
    def phase_volume_fractions(self):
        return np.array([i / self.num_voxels for i in self.phase_num_voxels])

    @property
    def phase_type_volume_fractions(self):
        return np.array([i / self.num_voxels for i in self.phase_type_num_voxels])

    @property
    def seeds_grid(self):
        return np.round(self.grid_size * self.seeds / self.size, decimals=0).astype(int)

    def get_voxel_IPF(self, IPF_dir=None, as_3D=False):
        if IPF_dir is None:
            IPF_dir = np.array([0, 0, 1])
        vox_oris = self.voxel_orientation
        dms_oris = Orientation(vox_oris.reshape((-1, 4)), family="cubic")
        IPF = dms_oris.IPF_color(IPF_dir)
        shape = list(vox_oris.shape[:-1]) + [3]
        vox_IPF = IPF.reshape(shape)

        if self.dimension == 2 and as_3D:
            return np.transpose(vox_IPF, axes=(1, 0, 2))[:, :, None, :]
        else:
            return vox_IPF

    def remove_interface(self, interface_name):
        """Remove an interface from the geometry. This will invalidate the geometry if
        the specified interface is referred by any phase-pairs."""

        idx = self.interface_names.index(interface_name)
        interface = self.interfaces.pop(idx)

        interface_map_tri = np.tril(-np.ones_like(self.interface_map)) + np.triu(
            self.interface_map
        )
        phase_pairs = np.array(np.where(interface_map_tri == idx))

        # set NaNs in interface map:
        self._interface_map[phase_pairs[0], phase_pairs[1]] = np.nan

        # realign indices in map that succeed the removed interface:
        self._interface_map[self._interface_map > idx] -= 1

        return interface, phase_pairs

    def get_grain_boundary_map(self, as_3D=False):
        voxel_GBs = np.ones_like(self.voxel_phase, dtype=int) * -10
        GBs = self.get_grain_boundaries().values()
        for idx, GB_i in enumerate(GBs):
            voxel_GBs[GB_i["voxel_indices"]] = idx

        if self.dimension == 3:
            return voxel_GBs
        elif as_3D:
            return voxel_GBs.T[:, :, None]

        return voxel_GBs

    def get_interface_energies_by_misorientation(self, misorientation_matrix=None):
        energies_theta = []
        if misorientation_matrix is None:
            misorientation_matrix = self.misorientation_matrix
        for interface_i in self.interfaces:
            pp_neighbours = []
            for pp in interface_i.phase_pairs:
                pp_is_neighbours = np.any(
                    np.all(pp[:, None] == self.neighbour_list, axis=0)
                )
                if pp_is_neighbours:
                    pp_neighbours.append(pp)
            pp_neighbours = np.array(pp_neighbours)
            if pp_neighbours.size:
                misoris = misorientation_matrix[pp_neighbours[:, 0], pp_neighbours[:, 1]]
                energies_theta.append(
                    {
                        "energy": interface_i.properties["energy"]["e0"],
                        "mobility": interface_i.properties["mobility"]["m0"],
                        "misorientation": misoris,
                        "phase_pairs": pp_neighbours,
                    }
                )
        return energies_theta

    def show_interface_energies_by_misorientation(
        self,
        interface_binning,
        misorientation_matrix=None,
        colour_by_bin=False,
        layout_args=None,
        show_mobility=True,
    ):
        energy_range = interface_binning.get("energy_range")
        mobility_range = interface_binning.get("mobility_range")

        if mobility_range is None:
            M_min, M_max = None, None
        else:
            M_min, M_max = mobility_range

        if energy_range is None:
            E_min, E_max = None, None
        else:
            E_min, E_max = energy_range

        fig = generate_interface_energies_plot(
            E_min=E_min,
            E_max=E_max,
            M_min=M_min,
            M_max=M_max,
            theta_max=interface_binning["theta_max"],
            n=interface_binning.get("n"),
            B=interface_binning.get("B"),
        )
        e_y = []
        m_y = []
        x = []
        color = []
        hover = []
        for bin_idx, bin_i in enumerate(
            self.get_interface_energies_by_misorientation(misorientation_matrix)
        ):
            for m_i_idx, m_i in enumerate(bin_i["misorientation"]):
                e_y.append(bin_i["energy"])
                m_y.append(bin_i["mobility"])
                x.append(m_i)
                color.append(bin_idx)
                hover.append(
                    f"({bin_i['phase_pairs'][m_i_idx, 0], bin_i['phase_pairs'][m_i_idx, 1]})"
                )

        fig.add_scatter(
            x=x,
            y=e_y,
            secondary_y=False,
            mode="markers",
            marker_color=color if colour_by_bin else "blue",
            marker_size=5,
            hovertext=hover,
            name="Model GB energies",
        )
        if show_mobility:
            fig.add_scatter(
                x=x,
                y=m_y,
                secondary_y=True,
                mode="markers",
                marker_color=color if colour_by_bin else "red",
                marker_size=5,
                hovertext=hover,
            )
        fig.layout.title = f"Interface properties at increment {self.increment}"
        fig.update_layout(layout_args or {})
        return fig

    def show_relative_misorientation(self, layout_args=None):
        phase_centroids = self.get_phase_voxel_centroids()
        idx = np.argmin(phase_centroids[:, 0])
        fig = px.imshow(
            np.rad2deg(
                quat_angle_between(
                    np.tile(
                        self.phase_orientation[idx], (self.phase_orientation.shape[0], 1)
                    ),
                    self.phase_orientation,
                )
            )[self.voxel_phase].T
        )
        fig.layout.update(layout_args or {})
        return fig
