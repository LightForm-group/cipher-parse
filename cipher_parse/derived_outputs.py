import numpy as np


def num_voxels_per_phase(input_data, phaseid):
    num_initial_phases = input_data["num_phases"]
    num_voxels_per_phase = np.zeros(num_initial_phases, dtype=int)
    uniq, counts = np.unique(phaseid.astype(int), return_counts=True)
    num_voxels_per_phase[uniq] = counts
    return num_voxels_per_phase
