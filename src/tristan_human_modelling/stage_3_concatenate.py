import os
import logging
from pathlib import Path


import pydmr
from miblab.pipe import run_stage

import utils.models_one_scan
import utils.models_two_scan

MODELS = {
    1: utils.models_one_scan,
    2: utils.models_two_scan,
}


def run(build, logfile):
    dir_output = os.path.join(build, 'human_modelling', 'stage_3_concatenate')
    dir_fits = os.path.join(build, 'human_modelling', 'stage_2_modelfit')

    for scans in [1, 2]:
        dir_pars = os.path.join(dir_fits, f"scans_{scans}", 'Pars')
        all_results_dmr = os.path.join(dir_output, f"scans_{scans}", 'all_results')
        all_results_wide = os.path.join(dir_output, f"scans_{scans}", 'all_results_wide.csv')
        all_sdev_wide = os.path.join(dir_output, f"scans_{scans}", 'all_results_sdev_wide.csv')
        
        subj_results = list(Path(dir_pars).glob("*.zip"))
        subj_results = [str(s) for s in subj_results]
        pydmr.concat(subj_results, all_results_dmr)
        pydmr.pars_to_wide(all_results_dmr, all_results_wide)
        pydmr.sdev_to_wide(all_results_dmr, all_sdev_wide)


if __name__ == '__main__':

    build = r"C:\Users\md1spsx\Documents\Data\tristan"
    run_stage(run, build, 'human_modelling', __file__)