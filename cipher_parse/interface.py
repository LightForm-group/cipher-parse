import copy
from typing import Dict, List, Optional, Tuple, Union

import numpy as np


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
