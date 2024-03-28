import copy
import json
import shutil
from subprocess import run, PIPE
from pathlib import Path
import re
import os
from textwrap import dedent

import numpy as np
import pyvista as pv
import pandas as pd
import plotly.express as px
import zarr

from cipher_parse.cipher_input import CIPHERInput, decompress_1D_array_string
from cipher_parse.geometry import CIPHERGeometry
from cipher_parse.utilities import (
    get_subset_indices,
    get_time_linear_subset_indices,
    update_plotly_figure_animation_slider_to_times,
)
from cipher_parse.derived_outputs import num_voxels_per_phase

DEFAULT_PARAVIEW_EXE = "pvbatch"
INC_DATA_NON_ARRAYS = (
    "increment",
    "time",
    "dimensions",
    "spacing",
    "number_VTI_cells",
    "number_VTI_points",
)
DERIVED_OUTPUTS_REQUIREMENTS = {
    "num_voxels_per_phase": ["phaseid"],
}
DERIVED_OUTPUTS_FUNCS = {
    "num_voxels_per_phase": num_voxels_per_phase,
}
STANDARD_OUTPUTS_TYPES = {
    "phaseid": int,
    "interfaceid": int,
    "matid": int,
}


def parse_cipher_stdout(path_or_string, is_string=False):
    warning_start = "Warning: "
    write_out = "writing output at time "

    warnings = []
    steps = []
    is_accepted = []
    time = []
    dt = []
    wlte = []
    wltea = []
    wlter = []
    outputs = {}  # keys file names; values times

    if is_string:
        lines = path_or_string.split("\n")
    else:
        with Path(path_or_string).open("rt") as fp:
            lines = fp.readlines()

    for ln_idx, ln in enumerate(lines):
        ln = ln.strip()
        if ln.startswith(warning_start):
            warnings.append(ln.split(warning_start)[1])
            continue

        step_search = re.search(r"\s+step\s+(\d+)\s+(.*)", ln)
        if step_search:
            groups = step_search.groups()
            step = int(groups[0])
            steps.append(step)

            step_dat = groups[1].split()

            is_accepted.append(bool(step_dat[0]))
            time.append(float(step_dat[1][2:].rstrip("+")))

            dt_pat = r"dt=(\d\.\d+e(-|\+)\d+)"
            dt_group = re.search(dt_pat, ln).groups()[0]
            dt.append(float(dt_group))

            wlte.append(float(step_dat[-5].lstrip("wlte=")))
            wltea.append(float(step_dat[-3]))
            wlter.append(float(step_dat[-1]))

        elif ln.startswith(write_out):
            ln_s = ln.split()
            outputs.update({ln_s[6]: float(ln_s[4])})
        else:
            continue

    out = {
        "warnings": warnings,
        "steps": np.array(steps),
        "is_accepted": np.array(is_accepted),
        "time": np.array(time),
        "dt": np.array(dt),
        "wlte": np.array(wlte),
        "wltea": np.array(wltea),
        "wlter": np.array(wlter),
        "outputs": outputs,
    }
    return out


def generate_VTI_files_from_VTU_files(
    sampling_dimensions,
    paraview_exe=DEFAULT_PARAVIEW_EXE,
):
    """Generate a 'ParaView-python' script for generating VTI files from VTU files and
    execute that script."""

    script_name = "vtu2vti.py"
    if len(sampling_dimensions) == 2:
        sampling_dimensions += [1]

    with Path(script_name).open("wt") as fp:
        fp.write(
            dedent(
                f"""
            import os

            from paraview.simple import *

            vtu_files = []
            for root, dirs, files in os.walk("."):
                for f in files:
                    if f.endswith(".vtu"):
                        vtu_files.append(f)
            
            for file_i_path in vtu_files:
                file_i_base_name = file_i_path.split(".")[0]
                vtu_data_i = XMLUnstructuredGridReader(
                    FileName=[os.getcwd() + os.path.sep + file_i_path]
                )
                resampleToImage1 = ResampleToImage(Input=vtu_data_i)
                resampleToImage1.SamplingDimensions = {sampling_dimensions!r}
                SetActiveSource(resampleToImage1)
                SaveData(file_i_base_name + ".vti", resampleToImage1)
        """
            )
        )

    proc = run(f"{paraview_exe} {script_name}", shell=True, stdout=PIPE, stderr=PIPE)
    stdout = proc.stdout.decode()
    stderr = proc.stderr.decode()
    if stdout:
        print(stdout)
    if stderr:
        print(stderr)


class CIPHEROutput:
    """Class to hold output information from a CIPHER simulation."""

    def __init__(
        self,
        directory,
        options,
        input_YAML_file_name,
        stdout_file_name,
        input_YAML_file_str,
        stdout_file_str,
        incremental_data,
        input_map_voxel_phase=None,
        input_map_phase_material=None,
        input_map_interface=None,
        quiet=False,
        cipher_input=None,
    ):
        default_options = {
            "paraview_exe": DEFAULT_PARAVIEW_EXE,
            "delete_VTIs": True,
            "delete_VTUs": False,
            "use_existing_VTIs": False,
            "num_VTU_files": None,
            "VTU_files_time_interval": None,
            "derive_outputs": None,
            "save_outputs": None,
        }

        self.directory = Path(directory)
        self.options = {**default_options, **options}
        self.input_YAML_file_name = input_YAML_file_name
        self.input_YAML_file_str = input_YAML_file_str
        self.input_map_voxel_phase = input_map_voxel_phase
        self.input_map_phase_material = input_map_phase_material
        self.input_map_interface = input_map_interface
        self.stdout_file_name = stdout_file_name
        self.stdout_file_str = stdout_file_str
        self.incremental_data = incremental_data
        self.quiet = quiet

        self._cipher_input = cipher_input or None
        self._cipher_stdout = None
        self._geometries = None  # assigned by set_geometries

        if (
            options.get("VTU_files_time_interval") is not None
            and options.get("num_VTU_files") is not None
        ):
            raise ValueError(
                "Specify at most one of 'num_VTU_files' and 'VTU_files_time_interval'."
            )

        for idx, i in enumerate(options["save_outputs"]):
            if i.get("number") is not None and i.get("time_interval") is not None:
                raise ValueError(
                    f"Specify at most one of 'number' and 'time_interval' for save "
                    f"output {idx}."
                )

    @classmethod
    def parse(
        cls,
        directory,
        options=None,
        input_YAML_file_name="cipher_input.yaml",
        stdout_file_name="stdout.log",
        get_voxel_phase=True,
        get_phase_material=True,
        get_interface=True,
    ):
        directory = Path(directory)

        yaml_path = directory / input_YAML_file_name
        with yaml_path.open("rt") as fp:
            input_YAML_file_str = "".join(fp.readlines())

        stdout_path = directory / stdout_file_name
        with stdout_path.open("rt") as fp:
            stdout_file_str = "".join(fp.readlines())

        (
            voxel_phase,
            phase_material,
            interface_map,
        ) = CIPHERInput.get_input_maps_from_files(
            inp_file_str=input_YAML_file_str,
            directory=directory,
            get_voxel_phase=get_voxel_phase,
            get_phase_material=get_phase_material,
            get_interface=get_interface,
        )

        obj = cls(
            directory=directory,
            options=options,
            input_YAML_file_name=input_YAML_file_name,
            input_YAML_file_str=input_YAML_file_str,
            input_map_voxel_phase=voxel_phase,
            input_map_phase_material=phase_material,
            input_map_interface=interface_map,
            stdout_file_name=stdout_file_name,
            stdout_file_str=stdout_file_str,
            incremental_data=None,
        )

        inc_data, outputs_keep_idx = obj.get_incremental_data()
        obj.incremental_data = inc_data
        obj.options["outputs_keep_idx"] = outputs_keep_idx

        return obj

    def _get_time_linear_subset_indices(self, time_interval):
        return get_time_linear_subset_indices(
            time_interval=time_interval,
            max_time=self.get_input_YAML_data()["solution_parameters"]["time"],
            times=np.array(list(self.cipher_stdout["outputs"].values())),
        )

    def get_incremental_data(self):
        """Generate temporary VTI files to parse requested cipher outputs on a uniform
        grid."""

        inp_dat = self.get_input_YAML_data()
        grid_size = inp_dat["grid_size"]

        if not self.options["use_existing_VTIs"]:
            generate_VTI_files_from_VTU_files(
                sampling_dimensions=grid_size,
                paraview_exe=self.options["paraview_exe"],
            )

        outfile_base = inp_dat["solution_parameters"]["outfile"]
        output_lookup = {
            i: f"{outfile_base} output.{idx}"
            for idx, i in enumerate(inp_dat["header"]["outputs"])
        }
        vtu_file_list = sorted(
            list(self.directory.glob(f"{outfile_base}_*.vtu")),
            key=lambda x: int(re.search(r"\d+", x.name).group()),
        )
        vti_file_list = sorted(
            list(self.directory.glob(f"{outfile_base}_*.vti")),
            key=lambda x: int(re.search(r"\d+", x.name).group()),
        )

        # Move all VTU files to a sub-directory:
        viz_dir = self.directory / "original_viz"

        vtu_orig_file_list = []
        if not viz_dir.is_dir():
            viz_dir.mkdir()
            for viz_file_i in vtu_file_list:
                dst_i = viz_dir.joinpath(viz_file_i.name).with_suffix(
                    ".viz" + viz_file_i.suffix
                )
                shutil.move(viz_file_i, dst_i)
                vtu_orig_file_list.append(dst_i)
        else:
            vtu_orig_file_list = sorted(
                list(viz_dir.glob("*")),
                key=lambda x: int(re.search(r"\d+", x.name).group()),
            )

        # Copy back to the root directory VTU files that we want to keep:
        if self.options["num_VTU_files"]:
            viz_files_keep_idx = get_subset_indices(
                len(vti_file_list),
                self.options["num_VTU_files"],
            )
        elif self.options["VTU_files_time_interval"]:
            viz_files_keep_idx = self._get_time_linear_subset_indices(
                time_interval=self.options["VTU_files_time_interval"]
            )
        else:
            viz_files_keep_idx = []

        for i in viz_files_keep_idx:
            viz_file_i = vtu_orig_file_list[i]
            dst_i = (
                self.directory.joinpath(viz_file_i.name)
                .with_suffix("")
                .with_suffix(".vtu")
            )
            shutil.copy(viz_file_i, dst_i)

        if self.options["delete_VTUs"]:
            print(f"Deleting original VTU files in directory: {viz_dir}")
            shutil.rmtree(viz_dir)

        # get which files to include for each output/derived output
        outputs_keep_idx = {}
        for save_out_i in self.options["save_outputs"]:
            if "number" in save_out_i:
                keep_idx = get_subset_indices(len(vti_file_list), save_out_i["number"])
            elif "time_interval" in save_out_i:
                keep_idx = self._get_time_linear_subset_indices(
                    time_interval=save_out_i["time_interval"]
                )
            else:
                keep_idx = list(range(len(vti_file_list)))
            outputs_keep_idx[save_out_i["name"]] = keep_idx

        incremental_data = []
        for file_i_idx, file_i in enumerate(vti_file_list):
            mesh = pv.get_reader(file_i).read()
            vtu_file_name = file_i.name.replace("vti", "vtu")
            inc_data_i = {
                "increment": int(re.search(r"\d+", file_i.name).group()),
                "time": self.cipher_stdout["outputs"][vtu_file_name],
                "dimensions": list(mesh.dimensions),
                "spacing": list(mesh.spacing),
                "number_VTI_cells": mesh.number_of_cells,
                "number_VTI_points": mesh.number_of_points,
            }

            standard_outputs = {}
            for name in output_lookup:
                arr_flat = mesh.get_array(output_lookup[name])
                arr = arr_flat.reshape(mesh.dimensions, order="F")
                if name in STANDARD_OUTPUTS_TYPES:
                    arr = arr.astype(STANDARD_OUTPUTS_TYPES[name])
                standard_outputs[name] = np.array(arr)  # convert from pyvista_ndarray

            derived_outputs = {}
            for derive_out_i in self.options["derive_outputs"]:
                name_i = derive_out_i["name"]
                func = DERIVED_OUTPUTS_FUNCS[name_i]
                func_args = {"input_data": inp_dat}
                func_args.update(
                    {i: standard_outputs[i] for i in DERIVED_OUTPUTS_REQUIREMENTS[name_i]}
                )
                derived_outputs[name_i] = func(**func_args)

            for out_name, keep_idx in outputs_keep_idx.items():
                if file_i_idx in keep_idx:
                    if out_name in DERIVED_OUTPUTS_REQUIREMENTS:
                        # a derived output:
                        inc_data_i[out_name] = derived_outputs[out_name]
                    else:
                        # a standard output:
                        inc_data_i[out_name] = standard_outputs[out_name]

            incremental_data.append(inc_data_i)

        if self.options["delete_VTIs"] and not self.options["use_existing_VTIs"]:
            for file_i in vti_file_list:
                print(f"Deleting temporary VTI file: {file_i}")
                os.remove(file_i)

        outputs_keep_idx["VTU_files"] = viz_files_keep_idx

        return incremental_data, outputs_keep_idx

    @property
    def cipher_input(self):
        if not self._cipher_input:
            self._cipher_input = CIPHERInput.from_input_YAML_str(
                file_str=self.input_YAML_file_str,
                quiet=self.quiet,
                input_map_voxel_phase=self.input_map_voxel_phase,
                input_map_phase_material=self.input_map_phase_material,
                input_map_interface=self.input_map_interface,
            )
        return self._cipher_input

    @property
    def cipher_stdout(self):
        if not self._cipher_stdout:
            if self.stdout_file_str:
                self._cipher_stdout = parse_cipher_stdout(
                    self.stdout_file_str,
                    is_string=True,
                )
            else:
                self._cipher_stdout = parse_cipher_stdout(
                    self.directory / self.stdout_file_name
                )
        return self._cipher_stdout

    def get_input_YAML_data(self, parse_interface_map=False):
        """Get some basic input details (using the YAML input file) without initialising
        the CIPHERInput object, which can take a while depending on the grid size."""
        dat = CIPHERInput.read_input_YAML_string(
            file_str=self.input_YAML_file_str,
            parse_interface_map=parse_interface_map,
        )
        if self.input_map_voxel_phase is not None:
            dat["voxel_phase"] = self.input_map_voxel_phase
            dat["unique_phase_IDs"] = np.unique(self.input_map_voxel_phase)

        if self.input_map_phase_material is not None:
            dat["phase_material"] = self.input_map_phase_material

        if self.input_map_interface is not None:
            dat["interface_map"] = self.input_map_interface

        return dat

    def to_JSON(self, keep_arrays=False):
        data = {
            "directory": str(self.directory),
            "options": self.options,
            "input_YAML_file_name": self.input_YAML_file_name,
            "input_YAML_file_str": self.input_YAML_file_str,
            "input_map_voxel_phase": self.input_map_voxel_phase,
            "input_map_phase_material": self.input_map_phase_material,
            "input_map_interface": self.input_map_interface,
            "stdout_file_name": self.stdout_file_name,
            "stdout_file_str": self.stdout_file_str,
            "incremental_data": self.incremental_data,
            "geometries": [i.to_JSON(keep_arrays) for i in self._geometries or []],
        }
        if not keep_arrays:
            for inc_idx, inc_i in enumerate(data["incremental_data"] or []):
                for key in inc_i:
                    if key not in INC_DATA_NON_ARRAYS:
                        as_list_val = np.copy(
                            data["incremental_data"][inc_idx][key]
                        ).tolist()
                        data["incremental_data"][inc_idx][key] = as_list_val

            if data["input_map_voxel_phase"] is not None:
                data["input_map_voxel_phase"] = np.copy(
                    data["input_map_voxel_phase"]
                ).tolist()

            if data["input_map_phase_material"] is not None:
                data["input_map_phase_material"] = np.copy(
                    data["input_map_phase_material"]
                ).tolist()

            if data["input_map_interface"] is not None:
                data["input_map_interface"] = np.copy(
                    data["input_map_interface"]
                ).tolist()

        return data

    @classmethod
    def from_JSON(cls, data, cipher_input=None, quiet=True):
        attrs = {
            "directory": data["directory"],
            "options": data["options"],
            "input_YAML_file_name": data["input_YAML_file_name"],
            "input_YAML_file_str": data["input_YAML_file_str"],
            "input_map_voxel_phase": data.get("input_map_voxel_phase"),
            "input_map_phase_material": data.get("input_map_phase_material"),
            "input_map_interface": data.get("input_map_interface"),
            "stdout_file_name": data["stdout_file_name"],
            "stdout_file_str": data["stdout_file_str"],
            "incremental_data": data["incremental_data"],
        }

        if attrs["input_map_voxel_phase"]:
            attrs["input_map_voxel_phase"] = np.array(attrs["input_map_voxel_phase"])

        if attrs["input_map_phase_material"]:
            attrs["input_map_phase_material"] = np.array(
                attrs["input_map_phase_material"]
            )

        if attrs["input_map_interface"]:
            attrs["input_map_interface"] = np.array(attrs["input_map_interface"])

        for inc_idx, inc_i in enumerate(attrs["incremental_data"] or []):
            for key, val in inc_i.items():
                if key not in INC_DATA_NON_ARRAYS and not isinstance(val, np.ndarray):
                    as_arr_val = np.array(attrs["incremental_data"][inc_idx][key])
                    attrs["incremental_data"][inc_idx][key] = as_arr_val

        obj = cls(**attrs, cipher_input=cipher_input, quiet=quiet)
        geoms = [
            CIPHERGeometry.from_JSON(i, quiet=quiet) for i in data.get("geometries", [])
        ]
        obj._geometries = geoms or None

        return obj

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

    def to_zarr(self, path):
        """Save to a persistent zarr store.

        This does not yet save `geometries`.

        """
        out_group = zarr.open_group(store=path)
        out_group.attrs.put(
            {
                "directory": str(self.directory),
                "options": self.options,
                "input_YAML_file_name": self.input_YAML_file_name,
                "input_map_voxel_phase": self.input_map_voxel_phase,
                "input_map_phase_material": self.input_map_phase_material,
                "input_map_interface": self.input_map_interface,
                "stdout_file_name": self.stdout_file_name,
            }
        )
        out_group.create_dataset(
            name="stdout_file_str",
            data=self.stdout_file_str.splitlines(),
        )
        out_group.create_dataset(
            name="input_YAML_file_str",
            data=self.input_YAML_file_str.splitlines(),
        )
        inc_dat_group = out_group.create_group("incremental_data", overwrite=True)
        for idx, inc_dat_i in enumerate(self.incremental_data):
            inc_dat_i_group = inc_dat_group.create_group(f"{idx}")
            inc_dat_i_group.attrs.put({k: inc_dat_i[k] for k in INC_DATA_NON_ARRAYS})
            for k in inc_dat_i:
                if k not in INC_DATA_NON_ARRAYS:
                    inc_dat_i_group.create_dataset(name=k, data=inc_dat_i[k])

        return out_group

    @classmethod
    def from_zarr(cls, path, cipher_input=None, quiet=True):
        """Load from a persistent zarr store.

        This does not yet load `geometries`.

        """
        group = zarr.open_group(store=path)
        attrs = group.attrs.asdict()
        kwargs = {
            "directory": attrs["directory"],
            "options": attrs["options"],
            "input_YAML_file_name": attrs["input_YAML_file_name"],
            "input_map_voxel_phase": attrs["input_map_voxel_phase"],
            "input_map_phase_material": attrs["input_map_phase_material"],
            "input_map_interface": attrs["input_map_interface"],
            "stdout_file_name": attrs["stdout_file_name"],
            "stdout_file_str": "\n".join(group.get("stdout_file_str")[:]),
            "input_YAML_file_str": "\n".join(group.get("input_YAML_file_str")[:]),
        }
        inc_data = []
        for inc_dat_i_group in group.get("incremental_data").values():
            inc_dat_i_group_attrs = inc_dat_i_group.attrs.asdict()
            inc_data_i = {k: inc_dat_i_group_attrs[k] for k in INC_DATA_NON_ARRAYS}
            for name, dataset in inc_dat_i_group.items():
                inc_data_i[name] = dataset[:]
            inc_data.append(inc_data_i)
        kwargs["incremental_data"] = inc_data

        obj = cls(**kwargs, cipher_input=cipher_input, quiet=quiet)
        return obj

    @classmethod
    def compare_phase_size_dist_evolution(
        cls,
        cipher_outputs,
        bin_size,
        use_phaseid=False,
        as_probability=False,
        max_increments=20,
        labels=None,
        row_labels=None,
        col_labels=None,
        label_order=None,
        row_label_name=None,
        col_label_name=None,
        label_name=None,
        layout_args=None,
    ):
        if len(cipher_outputs) > 1 and not (row_labels or col_labels or labels):
            raise TypeError(
                "Multiple cipher outputs but not labels/row_labels/col_labels "
                "specified."
            )

        if labels is not None:
            if len(labels) != len(cipher_outputs):
                raise TypeError(
                    "Length of `labels` must equal length of `cipher_outputs`."
                )
        elif not (row_labels or col_labels):
            labels = list(range(len(cipher_outputs)))

        if row_labels is not None:
            if len(row_labels) != len(cipher_outputs):
                raise TypeError(
                    "Length of `row_labels` must equal length of `cipher_outputs."
                )
        if col_labels is not None:
            if len(col_labels) != len(cipher_outputs):
                raise TypeError(
                    "Length of `col_labels` must equal length of `cipher_outputs."
                )

        label_name = label_name or "label"
        row_label_name = row_label_name or "row_label"
        col_label_name = col_label_name or "col_label"

        df_hist_all = pd.DataFrame()
        max_phase_size_all = 0
        max_prob_all = 0
        max_counts_all = 0
        for idx, out_i in enumerate(cipher_outputs):
            (
                df_hist_i,
                max_counts_i,
                max_prob_i,
                _,
                max_phase_size_i,
                available_inc_idx,
            ) = cls._prepare_phase_size_dist_evolution_dataframe(
                out_i,
                use_phaseid=use_phaseid,
                as_probability=as_probability,
                bin_size=bin_size,
                max_increments=max_increments,
            )

            if labels:
                df_hist_i[label_name] = str(labels[idx])
            if col_labels:
                df_hist_i[col_label_name] = str(col_labels[idx])
            if row_labels:
                df_hist_i[row_label_name] = str(row_labels[idx])

            df_hist_all = pd.concat([df_hist_all, df_hist_i])

            if max_phase_size_i > max_phase_size_all:
                max_phase_size_all = max_phase_size_i
            if max_prob_i > max_prob_all:
                max_prob_all = max_prob_i
            if max_counts_i > max_counts_all:
                max_counts_all = max_counts_i

        com_args = {}
        if row_labels:
            com_args["facet_row"] = row_label_name
        if col_labels:
            com_args["facet_col"] = col_label_name

        if label_order:
            com_args["category_orders"] = label_order

        if as_probability:
            if labels:
                com_args["color"] = label_name
            fig = px.bar(
                df_hist_all,
                x="bins",
                y="probability",
                labels={"x": "phase_size", "y": "probability"},
                animation_frame="evo_idx",
                barmode="overlay",
                **com_args,
            )
            y_max_lim = max_prob_all
        else:
            fig = px.bar(
                df_hist_all,
                x="bin",
                y="count",
                color="initial_bins",
                labels={"x": "phase_size", "y": "count"},
                animation_frame="evo_idx",
                barmode="overlay",
                **com_args,
            )
            y_max_lim = max_counts_all

        # turn off frame transitions:
        fig.layout.updatemenus[0].buttons[0].args[1]["transition"]["duration"] = 0

        fig.layout.update(
            {
                "xaxis": {
                    "range": [
                        -bin_size / 2,
                        np.round(max_phase_size_all * 1.1, decimals=6),
                    ],
                    "title": "phase size",
                },
                "yaxis": {"range": [0, y_max_lim]},
                "coloraxis": {
                    "colorbar": {"title": "Initial phase size"},
                    "colorscale": "viridis",
                },
                **(layout_args or {}),
            }
        )
        fig.update_traces(width=bin_size)
        fig.update_traces(marker_line={"width": 0})  # remove gap between stacked bars

        return fig

    @staticmethod
    def _prepare_phase_size_dist_evolution_dataframe(
        cipher_output,
        use_phaseid=False,
        as_probability=False,
        num_bins=None,
        bin_size=None,
        max_increments=20,
    ):
        input_yaml_dat = cipher_output.get_input_YAML_data()
        voxel_phase = input_yaml_dat["voxel_phase"]
        initial_phase_IDs = input_yaml_dat["unique_phase_IDs"]

        all_inc_data = cipher_output.incremental_data
        num_voxels_total = np.product(voxel_phase.shape)
        num_initial_phases = len(initial_phase_IDs)

        if use_phaseid:
            avail_inc_idx = cipher_output.options["outputs_keep_idx"]["phaseid"]
        else:
            avail_inc_idx = cipher_output.options["outputs_keep_idx"][
                "num_voxels_per_phase"
            ]

        subset_idx = get_subset_indices(len(avail_inc_idx), max_increments)
        avail_inc_idx = [avail_inc_idx[i] for i in subset_idx]

        num_incs = len(avail_inc_idx)

        if use_phaseid:
            num_voxels_per_phase = np.zeros((num_incs, num_initial_phases), dtype=int)
            for idx, inc_idx in enumerate(avail_inc_idx):
                inc_data = all_inc_data[inc_idx]
                phase_id = inc_data["phaseid"]
                uniq, counts = np.unique(phase_id, return_counts=True)
                num_voxels_per_phase[idx, uniq] = counts
        else:
            num_incs = len(avail_inc_idx)
            num_voxels_per_phase = np.vstack(
                [
                    all_inc_data[inc_idx]["num_voxels_per_phase"]
                    for inc_idx in avail_inc_idx
                ]
            )

        phase_size_normed = num_voxels_per_phase / num_voxels_total

        flattened_phase_size_normed = phase_size_normed.flatten()
        tiled_phase_ID = np.tile(np.arange(num_initial_phases), num_incs)
        repeated_incs = np.repeat(avail_inc_idx, num_initial_phases)

        # each row corresponds to a particular phase at a given increment:
        df = pd.DataFrame(
            {
                "phase_size": flattened_phase_size_normed,
                "phase_ID": tiled_phase_ID,
                "increment": repeated_incs,
            }
        )

        max_phase_size = df.phase_size.max()

        if num_bins is not None and bin_size is not None:
            raise TypeError(f"Specify exactly one of `num_bins` and `bin_size`.")
        elif num_bins is None and bin_size is None:
            num_bins = 50

        if bin_size is None:
            bin_size = max_phase_size / num_bins
        else:
            num_bins = int(max_phase_size / bin_size)

        bin_edges = np.linspace(0, max_phase_size + (bin_size / 2), num=num_bins + 1)
        bin_edges -= bin_size / 2  # so we have a bin centred on zero.

        df_hist = pd.DataFrame()
        initial_bins = None
        max_counts = 0
        max_prob = 0
        for evo_idx, inc_idx in enumerate(avail_inc_idx):
            df_inc_i = df[df["increment"] == inc_idx]
            counts, bins = np.histogram(df_inc_i.phase_size, bins=bin_edges)
            bin_centres = (np.array(bins) + (bin_size / 2))[:-1]
            prob_i = counts * bin_centres
            bin_indices_i = bins.searchsorted(
                df_inc_i.phase_size, "right"
            )  # bin index to which each phase belongs

            max_counts_i = np.max(counts)
            if max_counts_i > max_counts:
                max_counts = max_counts_i

            max_prob_i = np.max(prob_i)
            if max_prob_i > max_prob:
                max_prob = max_prob_i

            if as_probability:
                # each row corresponds to a single phase-size bin (for this increment):
                df_hist_i = pd.DataFrame(
                    {
                        "increment": [df_inc_i.increment.iat[0]] * num_bins,
                        "bins": bin_centres,
                        "probability": prob_i,
                    }
                )
            else:
                if inc_idx == 0:
                    initial_bins = bins[bin_indices_i - 1]

                # each row corresponds to a particular phase (for this increment):
                df_hist_i = pd.DataFrame(
                    {
                        "increment": df_inc_i.increment,
                        "initial_bins": initial_bins,
                        "phase_ID": df_inc_i.phase_ID,
                        "bin_index": bin_indices_i,
                        "bin": bins[bin_indices_i - 1],
                        "count": np.array([1] * len(bin_indices_i)),
                    }
                )

            df_hist_i["evo_idx"] = evo_idx

            df_hist = df_hist.append(df_hist_i)

        return df_hist, max_counts, max_prob, bin_size, max_phase_size, avail_inc_idx

    def show_phase_size_dist_evolution(
        self,
        use_phaseid=False,
        as_probability=False,
        num_bins=None,
        bin_size=None,
        max_increments=20,
        layout_args=None,
    ):
        """
        Parameters
        ----------
        use_phaseid : bool, optional
            If True, use the phaseid array to calculate the number of voxels per phase. If
            False, use the derived output `num_voxels_per_phase`.
        as_probability : bool, optional
            If True, the y-axis will be the probability of selecting a phase of a given
            size (binned number of voxels). If False, the y-axis will be simply the number
            of phases of a given size (binned number of voxels).
        layout_args : dict, optional
            Plotly layout options.
        """

        (
            df_hist,
            max_counts,
            max_prob,
            bin_size,
            max_phase_size,
            available_inc_idx,
        ) = self._prepare_phase_size_dist_evolution_dataframe(
            self,
            use_phaseid=use_phaseid,
            as_probability=as_probability,
            num_bins=num_bins,
            bin_size=bin_size,
            max_increments=max_increments,
        )

        if as_probability:
            fig = px.bar(
                df_hist,
                x="bins",
                y="probability",
                labels={"x": "phase_size", "y": "probability"},
                animation_frame="increment",
            )
            y_max_lim = max_prob
        else:
            fig = px.bar(
                df_hist,
                x="bin",
                y="count",
                color="initial_bins",
                labels={"x": "phase_size", "y": "count"},
                animation_frame="increment",
            )
            y_max_lim = max_counts

        # turn off frame transitions:
        fig.layout.updatemenus[0].buttons[0].args[1]["transition"]["duration"] = 0

        fig.layout.update(
            {
                "xaxis": {
                    "range": [
                        -bin_size / 2,
                        np.round(max_phase_size * 1.1, decimals=6),
                    ],
                    "title": "phase size",
                },
                "yaxis": {"range": [0, y_max_lim]},
                "coloraxis": {
                    "colorbar": {"title": "Initial phase size"},
                    "colorscale": "viridis",
                },
                **(layout_args or {}),
            }
        )
        fig.update_traces(width=bin_size)
        fig.update_traces(marker_line={"width": 0})  # remove gap between stacked bars

        times = [
            i["time"]
            for idx, i in enumerate(self.incremental_data)
            if idx in available_inc_idx
        ]
        update_plotly_figure_animation_slider_to_times(fig, times)

        return fig

    def show_misorientation_dist_evolution(
        self,
        num_bins=None,
        bin_size=None,
        layout_args=None,
    ):
        all_misori_vox = []
        inc_dat_indices = []
        incs = []
        times = []
        max_misori = 0
        for geom in self.geometries:
            misori_voxels = geom.voxel_map.get_interface_idx(
                self.cipher_input.geometry.misorientation_matrix
            ).flatten()
            misori_voxels = misori_voxels[misori_voxels != -1]
            all_misori_vox.append(misori_voxels)
            inc_dat_indices.append(geom.incremental_data_idx)
            incs.append(geom.increment)
            times.append(geom.time)
            max_misori_i = misori_voxels.max()
            if max_misori_i > max_misori:
                max_misori = max_misori_i

        if num_bins is not None and bin_size is not None:
            raise TypeError(f"Specify exactly one of `num_bins` and `bin_size`.")
        elif num_bins is None and bin_size is None:
            num_bins = 50

        if bin_size is None:
            bin_size = max_misori / num_bins
        else:
            num_bins = int(max_misori / bin_size)

        bin_edges = np.linspace(0, max_misori + (bin_size / 2), num=num_bins + 1)
        bin_edges -= bin_size / 2  # so we have a bin centred on zero.

        df_hist = pd.DataFrame()
        max_counts = 0
        for idx, voxels in enumerate(all_misori_vox):
            counts, bins = np.histogram(voxels, bins=bin_edges)
            max_counts_i = np.max(counts)
            if max_counts_i > max_counts:
                max_counts = max_counts_i
            bin_centres = (np.array(bins) + (bin_size / 2))[:-1]
            df_hist_i = pd.DataFrame(
                {
                    "incremental_data_idx": np.repeat(inc_dat_indices[idx], counts.size),
                    "increment": np.repeat(incs[idx], counts.size),
                    "time": np.repeat(times[idx], counts.size),
                    "misorientation": bin_centres,
                    "count": counts,
                }
            )
            df_hist = df_hist.append(df_hist_i)

        fig = px.bar(df_hist, x="misorientation", y="count", animation_frame="time")

        # turn off frame transitions:
        fig.layout.updatemenus[0].buttons[0].args[1]["transition"]["duration"] = 0

        fig.layout.update(
            {
                "xaxis": {
                    "range": [
                        -bin_size / 2,
                        np.round(max_misori * 1.1, decimals=6),
                    ],
                    "title": "Misorientation /degrees",
                },
                "yaxis": {"range": [0, max_counts], "title": "Num. voxels"},
                **(layout_args or {}),
            }
        )
        fig.update_traces(width=bin_size)
        fig.update_traces(marker_line={"width": 0})  # remove gap between stacked bars

        return fig

    def get_geometry(self, inc_data_index):
        start_geom = self.cipher_input.geometry
        inc_dat = self.incremental_data[inc_data_index]
        voxel_phase = inc_dat["phaseid"]
        if start_geom.dimension == 2:
            voxel_phase = voxel_phase[:, :, 0]

        geom = CIPHERGeometry(
            materials=start_geom.materials,
            interfaces=start_geom.interfaces,
            size=start_geom.size,
            voxel_phase=voxel_phase,
            allow_missing_phases=True,
            quiet=True,
            time=inc_dat["time"],
            increment=inc_dat["increment"],
            incremental_data_idx=inc_data_index,
        )
        return geom

    def get_all_geometries(self, include_initial=True):
        """A generator function to provide all available `CIPHERGeometry` objects."""

        if include_initial:
            geom_0 = self.cipher_input.geometry
            if geom_0.time is None:
                geom_0.time = 0
                geom_0.increment = 0
            yield geom_0

        if self._geometries is not None:
            for i in self._geometries:
                yield i
        else:
            for idx, inc_dat in enumerate(self.incremental_data):
                if "phaseid" in inc_dat:
                    yield self.get_geometry(idx)

    def set_all_geometries(self):
        if self._geometries is not None:
            raise ValueError("Geometries are already set.")
        self._geometries = [i for i in self.get_all_geometries(include_initial=True)]

    @property
    def geometries(self):
        if self._geometries is not None:
            return self._geometries
        else:
            raise ValueError("Run `set_all_geometries` first.")

    def show_slice_evolution(
        self,
        slice_index=0,
        normal_dir="z",
        data_label="phase",
        include=None,
        misorientation_matrix=None,
        layout_args=None,
        discrete_colours=None,
        step_size=None,
        **kwargs,
    ):
        """
        Parameters
        ----------
        discrete_colours : dict, optional
            If specified, a dict that maps slice data values to RGB three-tuples.
        """
        slices = []
        times = []
        step_size = step_size or 1
        for geom in self.geometries[::step_size]:
            slices.append(
                geom.get_slice(
                    slice_index, normal_dir, data_label, include, misorientation_matrix
                )[None]
            )
            times.append(geom.time)

        slices = np.concatenate(slices)
        if discrete_colours:
            slice_RGB = np.tile(slices[..., None], (1, 1, 1, 3)).astype(float)
            for k, v in discrete_colours.items():
                slice_RGB[np.where(np.all(slice_RGB == k, axis=3))] = v
            slices = slice_RGB
            fig = px.imshow(
                img=slices,
                animation_frame=0,
                **kwargs,
            )
        else:
            min_val = np.min(slices)
            max_val = np.max(slices)
            fig = px.imshow(
                img=slices,
                animation_frame=0,
                color_continuous_scale="viridis",
                zmin=min_val,
                zmax=max_val,
                labels={"color": data_label},
                **kwargs,
            )

        fig.update_layout(layout_args or {})
        update_plotly_figure_animation_slider_to_times(fig, times)
        return fig

    def get_average_radius_evolution(self, exclude=None):
        """Get an evolution proportional to the average radius across all phases."""

        inp_dat = self.get_input_YAML_data()
        num_phases = inp_dat["num_phases"]
        dimension = len(inp_dat["grid_size"])
        exclude = exclude or []

        num_voxel_inc_idx = []
        times = []
        for idx, i in enumerate(self.incremental_data):
            if "num_voxels_per_phase" in i:
                num_voxel_inc_idx.append(idx)
                times.append(i["time"])

        all_phases_num_voxels = np.ones((num_phases, len(num_voxel_inc_idx))) * np.nan
        for idx, inc_i_idx in enumerate(num_voxel_inc_idx):
            num_voxels = np.copy(self.incremental_data[inc_i_idx]["num_voxels_per_phase"])
            if exclude:
                num_voxels[exclude] = 0
            all_phases_num_voxels[:, idx] = num_voxels

        power = 1 / dimension
        all_phases_prop_radius = np.power(all_phases_num_voxels, power)
        all_phases_prop_radius[np.isclose(all_phases_prop_radius, 0)] = np.nan
        prop_avg_radius = np.nanmean(all_phases_prop_radius, axis=0)

        return prop_avg_radius, np.array(times)
