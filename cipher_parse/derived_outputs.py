import numpy as np


def num_voxels_per_phase(cipher_input, phaseid):
    voxel_phase = cipher_input.geometry.voxel_phase
    num_initial_phases = len(np.unique(voxel_phase))
    num_voxels_per_phase = np.zeros(num_initial_phases, dtype=int)
    uniq, counts = np.unique(phaseid.astype(int), return_counts=True)
    num_voxels_per_phase[uniq] = counts
    return num_voxels_per_phase
