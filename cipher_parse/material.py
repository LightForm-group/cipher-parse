import numpy as np

from cipher_parse.errors import (
    MaterialPhaseTypeFractionError,
    MaterialPhaseTypeLabelError,
    MaterialPhaseTypePhasesMissingError,
)


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

    @property
    def num_phases(self):
        return self.phases.size

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
