import json
from importlib import resources
import math
from pathlib import Path
from functools import reduce

import numpy as np
from scipy.spatial import Voronoi, Delaunay
from plotly import graph_objects
from plotly.colors import qualitative
from plotly.subplots import make_subplots
from vecmaths.geometry import get_box_xyz
from ipywidgets import interact


from cipher_parse.quats import axang2quat


def euclidean_distance_matrix(a, b):
    return np.linalg.norm(a[:, None, :] - b[None, :, :], axis=-1)


def in_hull(points, hull_points):
    """Test if a set of points are within a convex hull.

    Parameters
    ----------
    points : ndarray of shape (N, 3)
        Row vectors of points to test.
    hull_points : ndarray of shape (M, 3)
        Row vectors of vertices that are used to compute a convex hull

    """

    hull = Delaunay(hull_points)
    is_in = hull.find_simplex(points) >= 0

    return is_in


def get_coordinate_grid(size, grid_size):
    """Get the coordinates of the element centres of a uniform grid."""

    grid_size = np.array(grid_size)
    size = np.array(size)

    grid = np.meshgrid(*[np.arange(i) for i in grid_size])
    grid = np.moveaxis(np.array(grid), 0, -1)  # shape (*grid_size, dimension)

    element_size = (size / grid_size).reshape(1, 1, -1)

    coords = grid * size.reshape(1, 1, -1) / grid_size.reshape(1, 1, -1)
    coords += element_size / 2

    return coords, element_size


class DiscreteVoronoi:
    def __init__(self, seeds, grid_size, size=None, periodic=False, use_scipy=False):
        """
        Parameters
        ----------
        seeds : list or ndarray of shape (N, 2) or (N, 3)
            Row vectors of seed positions in 2D or 3D
        grid_size : list or ndarray of length 2 or 3
        size : list or ndarray of length 2 or 3, optional
            If not specified, a unit square/box is used.
        periodic : bool, optional
            Should the seeds and box be considered periodic. By default, False.

        """

        seeds = np.asarray(seeds)
        grid_size = np.asarray(grid_size)
        dimension = grid_size.size

        if size is None:
            size = [1] * dimension
        size = np.asarray(size)

        self.seeds = seeds
        self.grid_size = grid_size
        self.dimension = dimension
        self.size = size
        self.periodic = periodic

        self.element_size, self.coords, self.voxel_assignment = self._assign_voxels(
            use_scipy
        )

        self.coords_flat = self.coords.reshape(-1, self.dimension)
        self.voxel_assignment_flat = self.voxel_assignment.reshape(-1)

    def _assign_voxels(self, use_scipy):
        """Assign voxels to their closest seed point.

        Returns
        -------
        tuple
            element_size : ndarray of shape (`dimension`,)
            coords: ndarray of shape (*grid_size, `dimension`)
                Cartesian coordinates of the element centres
            voxel_assignment : ndarray of shape `grid_size`
                The index of the closest seed for each voxel.

        """

        # Get coordinates of grid element centres:
        coords, element_size = get_coordinate_grid(self.size, self.grid_size)

        # Get Euclidean distance of each grid element to each seed point
        coords_flat = coords.reshape(-1, self.dimension)

        seed_data = self.seeds_periodic if self.periodic else self.seeds

        if use_scipy:

            # Perform Voronoi tessellation of (periodic) seed positions:
            scipy_vor = Voronoi(seed_data)

            # Find which voxels are associated with each seed point:
            voxel_assignment = np.zeros(self.grid_size, dtype=int)
            for seed_idx, region_idx in enumerate(scipy_vor.point_region):
                region_verts_idx = scipy_vor.regions[region_idx]
                if region_verts_idx:
                    verts = scipy_vor.vertices[region_verts_idx]
                    voxels_in_flat = in_hull(coords_flat, verts)
                    voxels_in = voxels_in_flat.reshape(coords.shape[:-1])
                    voxels_in_idx = np.where(voxels_in)
                    voxel_assignment[voxels_in_idx] = seed_idx
                    # print(f'voxel assignment now: \n{voxel_assignment}')
        else:

            # shape (num coords, num seeds)
            dist = euclidean_distance_matrix(coords_flat, seed_data)  # BIG!

            # Assign each grid element to the nearest seed point
            nearest_seed_idx = np.argmin(dist, axis=1)  # shape (num coords,)
            voxel_assignment = nearest_seed_idx.reshape(coords.shape[:-1])

        if self.periodic:
            # Map periodic seed indices back to base seed indices:
            voxel_assignment = self.seeds_periodic_mapping[voxel_assignment]

        return element_size.flatten(), coords, voxel_assignment

    @property
    def num_seeds(self):
        return self.seeds.shape[0]

    @property
    def seeds_periodic_mapping(self):
        "Map periodic seed indices back to base seed indices."
        return np.tile(np.arange(self.seeds.shape[0]), 3**self.dimension)

    @property
    def seeds_periodic(self):
        """Get seeds positions including neighbouring periodic images."""
        trans_facts_2D = np.array(
            [
                [0, 0],
                [1, 0],
                [-1, 0],
                [0, 1],
                [0, -1],
                [1, 1],
                [-1, 1],
                [-1, -1],
                [1, -1],
            ]
        )
        if self.dimension == 2:
            trans = self.size * trans_facts_2D
        else:
            trans = self.size * np.vstack(
                [
                    np.hstack([trans_facts_2D, -np.ones((trans_facts_2D.shape[0], 1))]),
                    np.hstack([trans_facts_2D, np.zeros((trans_facts_2D.shape[0], 1))]),
                    np.hstack([trans_facts_2D, np.ones((trans_facts_2D.shape[0], 1))]),
                ]
            )

        seeds_periodic = np.concatenate(trans[:, None] + self.seeds)
        return seeds_periodic

    def show(self, show_voxels=False, show_periodic_seeds=False):
        """Show the discrete Voronoi tessellation."""

        if self.dimension == 2:
            edge_vecs = np.diag(np.concatenate([self.size, [0]]))
            box_xyz = get_box_xyz(edge_vecs, faces=True)["face01a"][0][[0, 1]]
        else:
            edge_vecs = np.diag(self.size)
            box_xyz = get_box_xyz(edge_vecs)[0]

        region_boundaries_approx_idx = np.where(
            np.diff(self.voxel_assignment, axis=0) != 0
        )
        region_boundaries_approx = self.coords[region_boundaries_approx_idx]
        print(region_boundaries_approx.shape)

        seed_data = self.seeds_periodic if show_periodic_seeds else self.seeds

        if self.dimension == 2:
            data = [
                {
                    "x": box_xyz[0],
                    "y": box_xyz[1],
                    "mode": "lines",
                    "name": "Box",
                },
                {
                    "x": region_boundaries_approx[:, 0],
                    "y": region_boundaries_approx[:, 1],
                    "name": "Boundaries approx.",
                    "mode": "markers",
                },
                {
                    "x": seed_data[:, 0],
                    "y": seed_data[:, 1],
                    "mode": "markers",
                    "name": "Seeds",
                    "marker": {
                        "size": 10,
                    },
                },
            ]
            if show_voxels:
                data.append(
                    {
                        "type": "scattergl",
                        "x": self.coords_flat[:, 0],
                        "y": self.coords_flat[:, 1],
                        "mode": "markers",
                        "marker": {
                            "color": self.voxel_assignment_flat,
                            "colorscale": qualitative.D3,
                            "size": 4,
                        },
                        "name": "Voxels",
                    }
                )

            layout = {
                "xaxis": {
                    "scaleanchor": "y",
                    "showgrid": False,
                    # 'dtick': self.element_size[0],
                },
                "yaxis": {
                    "showgrid": False,
                    # 'dtick': self.element_size[1],
                },
            }
        else:
            data = [
                {
                    "type": "scatter3d",
                    "x": box_xyz[0],
                    "y": box_xyz[1],
                    "z": box_xyz[2],
                    "mode": "lines",
                    "name": "Box",
                },
                {
                    "type": "scatter3d",
                    "x": region_boundaries_approx[:, 0],
                    "y": region_boundaries_approx[:, 1],
                    "z": region_boundaries_approx[:, 2],
                    "name": "Boundaries approx.",
                    "mode": "markers",
                    "marker": {
                        "size": 2,
                    },
                },
                {
                    "type": "scatter3d",
                    "x": seed_data[:, 0],
                    "y": seed_data[:, 1],
                    "z": seed_data[:, 2],
                    "mode": "markers",
                    "name": "Seeds",
                    "marker": {
                        "size": 10,
                    },
                },
            ]
            if show_voxels:
                data.append(
                    {
                        "type": "scatter3d",
                        "x": self.coords_flat[:, 0],
                        "y": self.coords_flat[:, 1],
                        "z": self.coords_flat[:, 2],
                        "mode": "markers",
                        "marker": {
                            "color": self.voxel_assignment_flat,
                            "size": 4,
                        },
                        "name": "Voxels",
                    }
                )
            layout = {}

        fig = graph_objects.FigureWidget(
            data=data,
            layout=layout,
        )

        return fig


def write_MTEX_EBSD_file(coords, euler_angles, phases, filename):

    dimension = coords.shape[1]

    col_names = [
        "Phase",
        "x",
        "y",
        "z",
        "Euler1",
        "Euler2",
        "Euler3",
    ]
    if dimension == 2:
        col_names.pop(3)

    all_dat = np.hstack(
        [
            phases[:, None],
            coords,
            euler_angles,
        ]
    )
    header = ", ".join([f"{i}" for i in col_names])
    np.savetxt(
        fname=filename,
        header=header,
        X=all_dat,
        fmt=["%d"] + ["%20.17f"] * (len(col_names) - 1),
    )


def write_MTEX_JSON_file(data, filename):
    with Path(filename).open("w") as fh:
        json.dump(data, fh)


def unjsonify_dict(dct):
    converted_to_list = dct.pop("_converted_to_list")
    for k, v in dct.items():
        if k in converted_to_list:
            v = np.array(v)
            dct[k] = v
    return dct


def jsonify_dict(dct):
    converted_to_list = []
    for k, v in dct.items():
        if isinstance(v, np.ndarray):
            v = v.tolist()
            converted_to_list.append(k)
        dct[k] = v
    dct["_converted_to_list"] = converted_to_list
    return dct


def get_by_path(root, path):
    """Get a nested dict or list item according to its "key path"

    Parameters
    ----------
    root : dict or list
        Can be arbitrarily nested.
    path : list of str
        The address of the item to get within the `root` structure.

    Returns
    -------
    sub_data : any

    """

    sub_data = root
    for key in path:
        sub_data = sub_data[key]

    return sub_data


def set_by_path(root, path, value):
    """Set a nested dict or list item according to its "key path"

    Parmaeters
    ----------
    root : dict or list
        Can be arbitrarily nested.
    path : list of str
        The address of the item to set within the `root` structure.
    value : any
        The value to set.

    """

    sub_data = root
    for key_idx, key in enumerate(path[:-1], start=1):
        try:
            sub_data = sub_data[key]
        except KeyError:
            set_by_path(sub_data, (key,), {})
            sub_data = sub_data[key]

    sub_data[path[-1]] = value


def read_shockley(theta, E_max, theta_max, degrees=True):
    """Misorientation-grain-boundary-energy relationship for low-angle GBs."""

    if degrees:
        theta = np.deg2rad(theta)
        theta_max = np.deg2rad(theta_max)

    zero_idx = np.isclose(theta, 0)

    A_0 = 1 + np.log(theta_max)
    E_0 = E_max / theta_max
    E = np.zeros_like(theta)
    E[~zero_idx] = E_0 * theta[~zero_idx] * (A_0 - np.log(theta[~zero_idx]))
    E[theta > theta_max] = np.max(E[theta < theta_max])

    return E


def grain_boundary_mobility(theta, M_max, theta_max, n=4, B=5, degrees=True):
    """Misorientation-grain-boundary-mobility relationship for low-angle GBs."""

    if degrees:
        theta = np.deg2rad(theta)
        theta_max = np.deg2rad(theta_max)

    M = M_max * (1 - np.exp(-B * (theta / theta_max) ** n))

    return M


def get_example_data_path_dream3D_2D():
    with resources.path(
        "cipher_parse.example_data.dream3d.2D", "synthetic_d3d.dream3d"
    ) as p:
        return p


def get_example_data_path_dream3D_3D():
    with resources.path(
        "cipher_parse.example_data.dream3d.3D", "synthetic_d3d.dream3d"
    ) as p:
        return p


def get_subset_indices(size, subset_size):
    """Get a list of N indices that index as uniformly as possible a sequence of a given
    size, with the constraint that the indices must include the initial and final elements.

    Parameters
    -----------
    size : int
    subset_size : int

    Returns
    -------
    list of int

    """
    if subset_size == 0:
        idx = []
    elif subset_size == 1 or (subset_size == 2 and size == 1):
        idx = [0]
    elif subset_size == 2:
        idx = [0, size - 1]
    elif subset_size >= size:
        idx = list(range(0, size))
    else:
        size_s = size - 1
        subset_size_s = subset_size - 1
        ratio = size_s / subset_size_s
        larger_group_size = math.ceil(ratio)
        smaller_group_size = math.floor(ratio)
        num_larger_groups = int(round(subset_size_s * (ratio % 1), ndigits=0))
        num_smaller_groups = subset_size_s - num_larger_groups

        if num_larger_groups >= num_smaller_groups:
            more_freq_group = (num_larger_groups, larger_group_size)
            less_freq_group = (num_smaller_groups, smaller_group_size)
        else:
            more_freq_group = (num_smaller_groups, smaller_group_size)
            less_freq_group = (num_larger_groups, larger_group_size)

        group_sizes = []
        num_A = more_freq_group[0]
        num_B = less_freq_group[0]

        if less_freq_group[0] > 0:

            num_ratio = int(more_freq_group[0] / less_freq_group[0])
            for i in range(subset_size_s):
                if num_A == num_B == 0:
                    break

                sub_i = []
                if num_A >= num_ratio:
                    sub_i = [more_freq_group[1]] * num_ratio
                    num_A -= num_ratio
                elif num_A < num_ratio:
                    sub_i = [more_freq_group[1]] * num_A
                    num_A = 0

                if num_B >= 1:
                    sub_i += [less_freq_group[1]] * 1
                    num_B -= 1
                group_sizes.extend(sub_i)
        else:
            group_sizes = [smaller_group_size] * num_smaller_groups

        idx = [0]
        for i in group_sizes[:-1]:
            idx.append(i + idx[-1])
        idx += [size_s]

    return idx


def get_time_linear_subset_indices(time_interval, max_time, times):
    intervals = np.linspace(
        0,
        max_time,
        num=int(((max_time + time_interval) / time_interval)),
        endpoint=True,
    )
    return list(set(np.argmin(np.abs(times - intervals[:, None]), axis=1).tolist()))


def sample_from_orientations_gradient(phase_centroids, max_misorientation_deg):
    """Generate an orientation gradient where orientations rotate uniformly according to a
    phase_centroids coordinate."""

    coords = phase_centroids[:, 0]
    low_x, high_x = np.min(coords), np.max(coords)
    frac_x = (coords - low_x) / (high_x - low_x)
    rots_deg = np.linspace(0, max_misorientation_deg, num=coords.size, endpoint=True)
    rots = np.deg2rad(rots_deg)
    ori_range = np.array([axang2quat(axis=np.array([0, 0, 1]), angle=i) for i in rots])
    ori_idx = np.argsort(frac_x)
    return ori_range, ori_idx


def generate_interface_energies_plot(
    E_min=0,
    M_min=0,
    E_max=1,
    M_max=1,
    theta_max=50,
    n=4,
    B=5,
):

    degrees = True
    theta = np.linspace(0, theta_max)

    plot_energy = E_min is not None and E_max is not None
    if plot_energy:
        E = (
            read_shockley(
                theta,
                E_max=(E_max - E_min),
                theta_max=theta_max,
                degrees=degrees,
            )
            + E_min
        )
    plot_mobility = M_min is not None and M_max is not None
    if plot_mobility:
        M = (
            grain_boundary_mobility(
                theta,
                M_max=(M_max - M_min),
                theta_max=theta_max,
                n=n,
                B=B,
                degrees=degrees,
            )
            + M_min
        )

    fig = make_subplots(specs=[[{"secondary_y": True}]])
    if plot_energy:
        fig.add_scatter(
            x=theta, y=E, name="Read-Shockley", secondary_y=False, line_color="blue"
        )
    if plot_mobility:
        fig.add_scatter(x=theta, y=M, name="mobility", secondary_y=True, line_color="red")

    fig.update_xaxes(showline=True, linewidth=1, linecolor="black", mirror=True)
    fig.update_yaxes(showline=True, linewidth=1, linecolor="black", mirror=True)
    fig.update_layout(
        legend_y=0.1,
        legend_x=0.9,
        legend_xanchor="right",
        legend_orientation="v",
    )
    fig.update_layout(
        {
            "width": 600,
            "xaxis_title": "Misori. /degrees",
        }
    )
    fig.update_yaxes(
        tickformat=".1e",
        secondary_y=False,
        title="GB energy",
        color="blue" if plot_mobility else None,
    )
    if plot_mobility:
        fig.update_yaxes(
            tickformat=".1e",
            secondary_y=True,
            color="red",
        )
    return fig


def generate_energy_widget(
    E_min=0, M_min=0, E_max=1, M_max=1, theta_max=50, n=4, B=5, degrees=True
):

    fig = generate_interface_energies_plot(
        E_min=E_min,
        M_min=M_min,
        E_max=E_max,
        M_max=M_max,
        theta_max=theta_max,
        n=n,
        B=B,
    )
    theta = np.linspace(0, theta_max)

    @interact(
        n=(1.0, 50.0, 0.01),
        B=(0, 50.0, 0.01),
        E_min=(0, 1, 0.01),
        M_min=(0, 1, 0.01),
        E_max=(0, 1, 0.01),
        M_max=(0, 1, 0.01),
        theta_max=(0, 90, 0.01),
    )
    def update(
        E_min=E_min, M_min=M_min, E_max=E_max, M_max=M_max, theta_max=theta_max, n=n, B=B
    ):
        with fig.batch_update():
            E = (
                read_shockley(
                    theta, E_max=(E_max - E_min), theta_max=theta_max, degrees=degrees
                )
                + E_min
            )
            M = (
                grain_boundary_mobility(
                    theta,
                    M_max=(M_max - M_min),
                    n=n,
                    B=B,
                    theta_max=theta_max,
                    degrees=degrees,
                )
                + M_min
            )
            fig.data[0].y = E
            fig.data[1].y = M

    return fig


def get_array_edge_mask(arr):
    """Get a boolean mask array that is True at the edge elements of an array."""
    all_idx = np.indices(arr.shape)
    mask = np.zeros_like(arr)
    for dim_idx, dim_size in enumerate(arr.shape):
        dim_mask = np.logical_or(all_idx[dim_idx] == 0, all_idx[dim_idx] == dim_size - 1)
        mask = np.logical_or(mask, dim_mask)
    return mask


def update_plotly_figure_animation_slider_to_times(fig, times):

    ani_steps = list(fig.layout.sliders[0]["steps"])
    ani_steps_new = []
    for idx, i in enumerate(ani_steps):
        i["label"] = f"{round(times[idx]):_}"
        ani_steps_new.append(i)

    fig.update_layout(
        sliders=[
            {
                "currentvalue": {"prefix": "Time = ", "suffix": " s"},
                "steps": ani_steps_new,
            }
        ]
    )
