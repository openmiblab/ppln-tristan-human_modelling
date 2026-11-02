import os
import logging
from pathlib import Path


import pydmr
from miblab.pipe import run_stage


PIPELINE = 'human_modelling'


def run(build, logfile):
    dir_input = os.path.join(build, PIPELINE, 'stage_1_modelfit')
    dir_output = os.path.join(build, PIPELINE, 'stage_2_concatenate')
    
    for scans in [[1], [1,2]]:
        dir_scans = f"scans_{len(scans)}"
        dir_pars = os.path.join(dir_input, dir_scans, 'Pars')
        all_results_dmr = os.path.join(dir_output, dir_scans, 'all_results')
        all_results_wide = os.path.join(dir_output, dir_scans, 'all_results_wide.csv')
        all_sdev_wide = os.path.join(dir_output, dir_scans, 'all_results_sdev_wide.csv')

        # Into pydmr as option
        subj_results = list(Path(dir_pars).glob("*.dmr.zip"))
        subj_results = [str(s) for s in subj_results]

        pydmr.concat(subj_results, all_results_dmr)
        pydmr.pars_to_wide(all_results_dmr, all_results_wide)
        pydmr.sdev_to_wide(all_results_dmr, all_sdev_wide)


if __name__ == '__main__':

    build = r"C:\Users\md1spsx\Documents\Data\tristan"
    run_stage(run, build, PIPELINE, __file__)