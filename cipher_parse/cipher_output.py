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

from cipher_parse.cipher_input import CIPHERInput
from cipher_parse.utilities import get_evenly_spaced_subset
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


def parse_cipher_stdout(path):
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

    with Path(path).open("rt") as fp:
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
        paraview_exe,
        delete_VTIs,
        delete_VTUs,
        max_viz_files,
        derive_outputs,
        save_outputs,
    ):

        self.directory = Path(directory)
        self.options = options
        self.input_YAML_file_name = input_YAML_file_name
        self.input_YAML_file_str = input_YAML_file_str
        self.stdout_file_name = stdout_file_name
        self.stdout_file_str = stdout_file_str
        self.incremental_data = incremental_data
        self.paraview_exe = paraview_exe
        self.delete_VTIs = delete_VTIs
        self.delete_VTUs = delete_VTUs
        self.max_viz_files = max_viz_files
        self.derive_outputs = derive_outputs
        self.save_outputs = save_outputs

        self._cipher_input = None
        self._cipher_stdout = None

    @classmethod
    def parse(
        cls,
        directory,
        options=None,
        input_YAML_file_name="cipher_input.yaml",
        stdout_file_name="stdout.log",
        paraview_exe=DEFAULT_PARAVIEW_EXE,
        delete_VTIs=True,
        delete_VTUs=False,
        use_existing_VTIs=False,
        max_viz_files=None,
        derive_outputs=None,
        save_outputs=None,
    ):
        directory = Path(directory)

        yaml_path = directory / input_YAML_file_name
        with yaml_path.open("rt") as fp:
            input_YAML_file_str = "".join(fp.readlines())

        stdout_path = directory / stdout_file_name
        with stdout_path.open("rt") as fp:
            stdout_file_str = "".join(fp.readlines())

        obj = cls(
            directory=directory,
            options=options,
            input_YAML_file_name=input_YAML_file_name,
            input_YAML_file_str=input_YAML_file_str,
            stdout_file_name=stdout_file_name,
            stdout_file_str=stdout_file_str,
            paraview_exe=paraview_exe,
            incremental_data=None,
            delete_VTIs=delete_VTIs,
            delete_VTUs=delete_VTUs,
            max_viz_files=max_viz_files,
            derive_outputs=derive_outputs,
            save_outputs=save_outputs,
        )

        obj.incremental_data = obj.get_incremental_data(
            use_existing_VTIs=use_existing_VTIs,
            delete_VTIs=delete_VTIs,
            delete_VTUs=delete_VTUs,
            max_viz_files=max_viz_files,
            derive_outputs=derive_outputs,
            save_outputs=save_outputs,
        )

        if max_viz_files is not None:
            pass

            viz_dir = Path("original_viz")
            viz_dir.mkdir()

            vtu_file_list = sorted(
                list(obj.directory.glob("*.vtu")),
                key=lambda x: int(re.search(r"\d+", x.name).group()),
            )
            vti_file_list = sorted(
                list(obj.directory.glob("*.vti")),
                key=lambda x: int(re.search(r"\d+", x.name).group()),
            )

            vtu_orig_file_list = []
            for viz_file_i in vtu_file_list:
                dst_i = viz_dir.joinpath(viz_file_i.name).with_suffix(
                    ".viz" + viz_file_i.suffix
                )
                shutil.move(viz_file_i, dst_i)
                vtu_orig_file_list.append(dst_i)

            for viz_file_i in vti_file_list:
                dst_i = viz_dir.joinpath(viz_file_i.name).with_suffix(
                    ".viz" + viz_file_i.suffix
                )
                shutil.move(viz_file_i, dst_i)

            viz_files_idx = get_evenly_spaced_subset(vtu_file_list, max_viz_files)

            for i in viz_files_idx:
                viz_file_i = vtu_orig_file_list[i]
                dst_i = (
                    Path("").joinpath(viz_file_i.name).with_suffix("").with_suffix(".vtu")
                )
                shutil.copy(viz_file_i, dst_i)

        return obj

    def get_incremental_data(
        self,
        use_existing_VTIs,
        delete_VTIs,
        delete_VTUs,
        max_viz_files,
        derive_outputs,
        save_outputs,
    ):
        """Generate temporary VTI files to parse requested cipher outputs on a uniform
        grid."""

        cipher_input = self.cipher_input

        if not use_existing_VTIs:
            generate_VTI_files_from_VTU_files(
                sampling_dimensions=cipher_input.geometry.grid_size.tolist(),
                paraview_exe=self.paraview_exe,
            )

        output_lookup = {
            i: f"out output.{idx}" for idx, i in enumerate(self.cipher_input.outputs)
        }
        vtu_file_list = sorted(
            list(self.directory.glob("*.vtu")),
            key=lambda x: int(re.search(r"\d+", x.name).group()),
        )
        vti_file_list = sorted(
            list(self.directory.glob("*.vti")),
            key=lambda x: int(re.search(r"\d+", x.name).group()),
        )

        # Move all VTU files to a sub-directory:
        viz_dir = Path("original_viz")
        viz_dir.mkdir()
        vtu_orig_file_list = []
        for viz_file_i in vtu_file_list:
            dst_i = viz_dir.joinpath(viz_file_i.name).with_suffix(
                ".viz" + viz_file_i.suffix
            )
            shutil.move(viz_file_i, dst_i)
            vtu_orig_file_list.append(dst_i)

        # Copy back to the root directory VTU files that we want to keep:
        viz_files_keep_idx = get_evenly_spaced_subset(vti_file_list, max_viz_files)
        for i in viz_files_keep_idx:
            viz_file_i = vtu_orig_file_list[i]
            dst_i = Path("").joinpath(viz_file_i.name).with_suffix("").with_suffix(".vtu")
            shutil.copy(viz_file_i, dst_i)

        if delete_VTUs:
            print(f"Deleting temporary original VTU files in directory: {viz_dir}")
            # shutil.rmtree(viz_dir) # TEMP

        # get which files to include for each output/derived output
        outputs_keep_idx = {}
        for save_out_i in save_outputs:
            if "max_num" in save_out_i:
                keep_idx = get_evenly_spaced_subset(vti_file_list, save_out_i["max_num"])
            else:
                keep_idx = range(len(vti_file_list))
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
                arr_flat = mesh.get_array(output_lookup[out_name])
                arr = arr_flat.reshape(mesh.dimensions, order="F")
                standard_outputs[name] = arr

            derived_outputs = {}
            for derive_out_i in derive_outputs:
                name_i = derive_out_i["name"]
                func = DERIVED_OUTPUTS_FUNCS[name_i]
                func_args = {"cipher_input": cipher_input}
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

        if delete_VTIs and not use_existing_VTIs:
            for file_i in vti_file_list:
                print(f"Deleting temporary VTI file: {file_i}")
                os.remove(file_i)

        return incremental_data

    @property
    def cipher_input(self):
        if not self._cipher_input:
            self._cipher_input = CIPHERInput.from_input_YAML_str(self.input_YAML_file_str)
        return self._cipher_input

    @property
    def cipher_stdout(self):
        if not self._cipher_stdout:
            self._cipher_stdout = parse_cipher_stdout(
                self.directory / self.stdout_file_name
            )
        return self._cipher_stdout

    def to_JSON(self, keep_arrays=False):
        data = {
            "directory": str(self.directory),
            "options": self.options,  # TODO
            "input_YAML_file_name": self.input_YAML_file_name,
            "input_YAML_file_str": self.input_YAML_file_str,
            "stdout_file_name": self.stdout_file_name,
            "stdout_file_str": self.stdout_file_str,
            "incremental_data": self.incremental_data,
            "paraview_exe": self.paraview_exe,
        }
        if not keep_arrays:
            for inc_idx, inc_i in enumerate(data["incremental_data"] or []):
                for key in inc_i:
                    if key not in INC_DATA_NON_ARRAYS:
                        print(f"key is array-like: {key}; val type is {type(inc_i[key])}")
                        print(type(data["incremental_data"][inc_idx][key]))
                        as_list_val = data["incremental_data"][inc_idx][key].tolist()
                        data["incremental_data"][inc_idx][key] = as_list_val

        return data

    @classmethod
    def from_JSON(cls, data):

        attrs = {
            "directory": data["directory"],
            "options": data["options"],
            "input_YAML_file_name": data["input_YAML_file_name"],
            "input_YAML_file_str": data["input_YAML_file_str"],
            "stdout_file_name": data["stdout_file_name"],
            "stdout_file_str": data["stdout_file_str"],
            "incremental_data": data["incremental_data"],
            "paraview_exe": data["paraview_exe"],
        }

        for inc_idx, inc_i in enumerate(attrs["incremental_data"] or []):
            for key, val in inc_i.items():
                if key not in INC_DATA_NON_ARRAYS and not isinstance(val, np.ndarray):
                    as_arr_val = np.array(attrs["incremental_data"][inc_idx][key])
                    attrs["incremental_data"][inc_idx][key] = as_arr_val

        obj = cls(**attrs)
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

    def get_num_voxels_per_phase(self):
        pass

    def show_phase_size_dist_evolution(self):

        voxel_phase = self.cipher_input.geometry.voxel_phase
        all_inc_data = self.incremental_data

        initial_phase_IDs = np.unique(voxel_phase)
        num_initial_phases = len(initial_phase_IDs)
        num_voxels_total = np.product(voxel_phase.shape)
        num_incs = len(all_inc_data)
        num_voxels_per_phase = np.zeros((num_incs, num_initial_phases), dtype=int)

        for inc_idx, inc_data in enumerate(all_inc_data):
            phase_id = inc_data["phaseid"]
            uniq, counts = np.unique(phase_id.astype(int), return_counts=True)
            num_voxels_per_phase[inc_idx, uniq] = counts

        phase_size_normed = num_voxels_per_phase / num_voxels_total
        flattened_phase_size_normed = phase_size_normed.flatten()
        tiled_phase_ID = np.tile(np.arange(num_initial_phases), num_incs)
        repeated_incs = np.repeat(np.arange(num_incs), num_initial_phases)
        df = pd.DataFrame(
            {
                "phase_size": flattened_phase_size_normed,
                "phase_ID": tiled_phase_ID,
                "increment": repeated_incs,
            }
        )
        num_bins = 50
        bin_size = 4 / (num_bins * num_initial_phases)
        bin_edges = np.linspace(0, df.phase_size.max(), num=num_bins)

        df_hist = pd.DataFrame()
        initial_bins = None
        max_counts = 0
        for inc_idx, _ in enumerate(all_inc_data):
            df_inc_i = df[df["increment"] == inc_idx]
            counts, bins = np.histogram(df_inc_i.phase_size, bins=bin_edges)
            bin_indices_i = bins.searchsorted(
                df_inc_i.phase_size, "right"
            )  # bin index to which each phase belongs

            max_counts_i = np.max(counts)
            if max_counts_i > max_counts:
                max_counts = max_counts_i

            if inc_idx == 0:
                initial_bins = bins[bin_indices_i - 1]

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
            df_hist = df_hist.append(df_hist_i)

        fig = px.bar(
            df_hist,
            x="bin",
            y="count",
            color="initial_bins",
            labels={"x": "phase_size", "y": "count"},
            animation_frame="increment",
        )

        # turn off frame transitions:
        fig.layout.updatemenus[0].buttons[0].args[1]["transition"]["duration"] = 0

        fig.layout.update(
            {
                "xaxis": {
                    "range": [
                        -bin_size / 2,
                        np.round(np.max(flattened_phase_size_normed) * 1.1, decimals=2),
                    ],
                    "title": "phase size",
                },
                "yaxis": {"range": [0, max_counts]},
                "width": 600,
                "coloraxis": {
                    "colorbar": {"title": "Initial phase size"},
                    "colorscale": "viridis",
                },
            }
        )
        fig.update_traces(width=bin_size)

        return fig
