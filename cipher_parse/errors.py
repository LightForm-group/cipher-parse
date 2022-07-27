class MaterialPhaseTypeFractionError(Exception):
    pass


class MaterialPhaseTypeLabelError(Exception):
    pass


class GeometryMissingPhaseAssignmentError(Exception):
    pass


class GeometryVoxelPhaseError(Exception):
    pass


class GeometryUnassignedPhasePairInterfaceError(Exception):
    pass


class MaterialPhaseTypePhasesMissingError(Exception):
    pass


class GeometryDuplicateMaterialNameError(Exception):
    pass


class GeometryNonUnitTargetVolumeFractionError(Exception):
    pass


class GeometryExcessTargetVolumeFractionError(GeometryNonUnitTargetVolumeFractionError):
    pass
