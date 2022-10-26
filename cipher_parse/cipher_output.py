from pathlib import Path
import re

import numpy as np

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
    outputs = {} # keys file names; values times

    with Path(path).open('rt') as fp:
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
                
                wlte.append(float(step_dat[-5].lstrip('wlte=')))
                wltea.append(float(step_dat[-3]))
                wlter.append(float(step_dat[-1]))

            elif ln.startswith(write_out):
                ln_s = ln.split()
                outputs.update({
                    ln_s[6]: float(ln_s[4])
                })
            else:
                continue

    out = {
        'warnings': warnings,
        'steps': np.array(steps),
        'is_accepted': np.array(is_accepted),
        'time': np.array(time),
        'dt': np.array(dt),
        'wlte': np.array(wlte),
        'wltea': np.array(wltea),
        'wlter': np.array(wlter),
        'outputs': outputs,

    }
    return out
